//! Self-consistent field: restricted and unrestricted Hartree-Fock, restricted Kohn-Sham DFT.

use crate::basis::{Basis, BasisName};
use crate::eri::{build_jk, compute_eris, EriTensor};
use crate::grid::{molecular_grid, GridLevel, MolecularGrid};
use crate::linalg::{congruence, eigh, matmul, solve, transpose};
use crate::molecule::{Atom, Molecule};
use crate::numint::xc_matrix;
use crate::one_electron::{dipole_matrices, one_electron_matrices};
use crate::xc::Functional;
use crate::Error;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    HartreeFock,
    Dft(Functional),
}

impl Method {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s.to_ascii_lowercase().as_str() {
            "hf" => Ok(Method::HartreeFock),
            "lda" | "svwn" => Ok(Method::Dft(Functional::Lda)),
            "pbe" => Ok(Method::Dft(Functional::Pbe)),
            _ => Err(Error::new(format!("Unknown method: {s}"))),
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            Method::HartreeFock => "hf",
            Method::Dft(Functional::Lda) => "lda",
            Method::Dft(Functional::Pbe) => "pbe",
        }
    }

    /// Fraction of exact exchange.
    pub fn exact_exchange(self) -> f64 {
        match self {
            Method::HartreeFock => 1.0,
            Method::Dft(f) => f.exact_exchange(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct ScfOptions {
    pub basis: BasisName,
    pub method: Method,
    pub grid: GridLevel,
    pub max_iterations: usize,
    /// Energy change threshold (Hartree).
    pub conv_energy: f64,
    /// Max |FDS - SDF| threshold in the orthogonal basis.
    pub conv_error: f64,
    /// Store orbitals after every iteration (for playback).
    pub record_history: bool,
}

impl ScfOptions {
    pub fn new(basis: BasisName, method: Method) -> Self {
        ScfOptions {
            basis,
            method,
            grid: GridLevel::Medium,
            max_iterations: 100,
            conv_energy: 1e-8,
            conv_error: 1e-5,
            record_history: false,
        }
    }
}

/// Orbitals from an earlier, nearby geometry used as the starting guess.
#[derive(Clone, Debug)]
pub struct OrbitalGuess {
    pub basis: BasisName,
    pub nbf: usize,
    pub nmo: usize,
    pub ca: Vec<f64>,
    pub cb: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct ScfIteration {
    pub energy: f64,
    pub error: f64,
    pub ca: Vec<f64>,
    pub cb: Vec<f64>,
}

pub struct ScfResult {
    pub converged: bool,
    pub energy: f64,
    pub iterations: usize,
    pub unrestricted: bool,
    pub method: Method,
    pub basis: Basis,
    pub nbf: usize,
    pub nmo: usize,
    pub nalpha: usize,
    pub nbeta: usize,
    /// MO coefficients, nbf x nmo row-major (column k is orbital k).
    pub ca: Vec<f64>,
    pub cb: Vec<f64>,
    pub eps_a: Vec<f64>,
    pub eps_b: Vec<f64>,
    pub da: Vec<f64>,
    pub db: Vec<f64>,
    pub nuclear_repulsion: f64,
    pub exc: f64,
    pub dipole: [f64; 3],
    pub history: Vec<ScfIteration>,
    pub grid: Option<MolecularGrid>,
}

impl ScfResult {
    pub fn guess(&self) -> OrbitalGuess {
        OrbitalGuess {
            basis: self.basis.name,
            nbf: self.nbf,
            nmo: self.nmo,
            ca: self.ca.clone(),
            cb: self.cb.clone(),
        }
    }

    pub fn total_density(&self) -> Vec<f64> {
        self.da.iter().zip(&self.db).map(|(a, b)| a + b).collect()
    }
}

fn density_from(c: &[f64], n: usize, nmo: usize, nocc: usize) -> Vec<f64> {
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = 0.0;
            for k in 0..nocc {
                s += c[i * nmo + k] * c[j * nmo + k];
            }
            d[i * n + j] = s;
            d[j * n + i] = s;
        }
    }
    d
}

/// C (C^T S C)^(-1/2)
fn orthonormalize(c: &[f64], s: &[f64], n: usize, m: usize) -> Vec<f64> {
    let mm = congruence(c, s, n, m);
    let e = eigh(&mm, m);
    let mut inv = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut acc = 0.0;
            for k in 0..m {
                acc += e.vectors[i * m + k] * e.vectors[j * m + k] / e.values[k].max(1e-14).sqrt();
            }
            inv[i * m + j] = acc;
        }
    }
    matmul(c, &inv, n, m, m)
}

type AtomCache = Mutex<HashMap<(u32, BasisName), Vec<f64>>>;

/// Spherically averaged atomic density from fractional-occupation restricted HF on the free atom.
fn atomic_density(z: u32, basis_name: BasisName) -> Result<Vec<f64>, Error> {
    static CACHE: OnceLock<AtomCache> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(d) = cache.lock().unwrap().get(&(z, basis_name)) {
        return Ok(d.clone());
    }
    let atoms = [Atom { z, xyz: [0.0; 3] }];
    let basis = Basis::new(&atoms, basis_name)?;
    let n = basis.nbf;
    let one = one_electron_matrices(&basis, &atoms);
    let h: Vec<f64> = one.t.iter().zip(&one.v).map(|(t, v)| t + v).collect();
    let eri = compute_eris(&basis);
    let se = eigh(&one.s, n);
    let mut x = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            x[i * n + k] = se.vectors[i * n + k] / se.values[k].sqrt();
        }
    }
    let xt = transpose(&x, n, n);
    let nelec = z as f64;
    let density_from_fock = |f: &[f64]| -> Vec<f64> {
        let e = eigh(&matmul(&xt, &matmul(f, &x, n, n, n), n, n, n), n);
        let c = matmul(&x, &e.vectors, n, n, n);
        let mut occ = vec![0.0; n];
        let mut left = nelec;
        let mut k = 0;
        while k < n && left > 1e-12 {
            let mut m = k + 1;
            while m < n && (e.values[m] - e.values[k]).abs() < 1e-4 {
                m += 1;
            }
            let put = left.min(2.0 * (m - k) as f64);
            for o in occ.iter_mut().take(m).skip(k) {
                *o = put / (m - k) as f64;
            }
            left -= put;
            k = m;
        }
        let mut d = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for (q, &o) in occ.iter().enumerate() {
                    if o > 0.0 {
                        s += o * c[i * n + q] * c[j * n + q];
                    }
                }
                d[i * n + j] = s;
            }
        }
        d
    };
    let mut d = density_from_fock(&h);
    for _ in 0..80 {
        let (j, k) = build_jk(&eri, &d, &[&d]);
        let f: Vec<f64> = (0..n * n).map(|q| h[q] + j[q] - 0.5 * k[0][q]).collect();
        let dn = density_from_fock(&f);
        let mut diff: f64 = 0.0;
        for q in 0..n * n {
            diff = diff.max((dn[q] - d[q]).abs());
            d[q] = 0.5 * d[q] + 0.5 * dn[q];
        }
        if diff < 1e-6 {
            break;
        }
    }
    cache.lock().unwrap().insert((z, basis_name), d.clone());
    Ok(d)
}

/// Superposition of atomic densities (total density).
fn sad_density(mol: &Molecule, basis: &Basis) -> Result<Vec<f64>, Error> {
    let n = basis.nbf;
    let mut d = vec![0.0; n * n];
    let mut offset = 0;
    for atom in &mol.atoms {
        let da = atomic_density(atom.z, basis.name)?;
        let m = (da.len() as f64).sqrt().round() as usize;
        for i in 0..m {
            for j in 0..m {
                d[(offset + i) * n + offset + j] = da[i * m + j];
            }
        }
        offset += m;
    }
    let ntot: i64 = mol.atoms.iter().map(|a| a.z as i64).sum();
    let nelec = mol.electron_count();
    if nelec != ntot {
        let f = nelec as f64 / ntot as f64;
        d.iter_mut().for_each(|v| *v *= f);
    }
    Ok(d)
}

struct DiisEntry {
    f: Vec<Vec<f64>>,
    e: Vec<f64>,
}

fn diis_extrapolate(entries: &[DiisEntry]) -> Vec<Vec<f64>> {
    let m = entries.len();
    let mut b = vec![0.0; (m + 1) * (m + 1)];
    for i in 0..m {
        for j in 0..=i {
            let s: f64 = entries[i]
                .e
                .iter()
                .zip(&entries[j].e)
                .map(|(x, y)| x * y)
                .sum();
            b[i * (m + 1) + j] = s;
            b[j * (m + 1) + i] = s;
        }
        b[i * (m + 1) + m] = -1.0;
        b[m * (m + 1) + i] = -1.0;
    }
    let mut rhs = vec![0.0; m + 1];
    rhs[m] = -1.0;
    match solve(&b, &rhs, m + 1) {
        Some(c) => (0..entries[0].f.len())
            .map(|s| {
                let size = entries[0].f[s].len();
                let mut f = vec![0.0; size];
                for (i, en) in entries.iter().enumerate() {
                    for k in 0..size {
                        f[k] += c[i] * en.f[s][k];
                    }
                }
                f
            })
            .collect(),
        None => entries[m - 1].f.clone(),
    }
}

/// Run an SCF calculation. `guess` (orbitals from a nearby geometry) speeds up convergence.
pub fn run_scf(
    mol: &Molecule,
    opts: &ScfOptions,
    guess: Option<&OrbitalGuess>,
) -> Result<ScfResult, Error> {
    let nelec = mol.electron_count();
    let mult = mol.multiplicity as i64;
    if mol.atoms.is_empty() {
        return Err(Error::new("The molecule has no atoms."));
    }
    if nelec <= 0 {
        return Err(Error::new("The molecule has no electrons."));
    }
    if mult < 1 || (nelec + mult - 1) % 2 != 0 || mult - 1 > nelec {
        return Err(Error::new(format!(
            "{nelec} electrons cannot have multiplicity {mult}."
        )));
    }
    let nalpha = ((nelec + mult - 1) / 2) as usize;
    let nbeta = nelec as usize - nalpha;
    let unrestricted = mult != 1;
    let xc = match opts.method {
        Method::HartreeFock => None,
        Method::Dft(f) => Some(f),
    };
    if xc.is_some() && unrestricted {
        return Err(Error::new(
            "Open-shell DFT is not supported yet; use Hartree-Fock (UHF).",
        ));
    }
    let kx = opts.method.exact_exchange();

    let basis = Basis::new(&mol.atoms, opts.basis)?;
    let n = basis.nbf;
    let one = one_electron_matrices(&basis, &mol.atoms);
    let h: Vec<f64> = one.t.iter().zip(&one.v).map(|(t, v)| t + v).collect();
    let s = one.s;
    let eri: EriTensor = compute_eris(&basis);
    let grid = xc.map(|_| molecular_grid(&mol.atoms, opts.grid));
    let enuc = mol.nuclear_repulsion();

    // Canonical orthogonalization X = U s^(-1/2), dropping near-linear dependencies.
    let se = eigh(&s, n);
    let keep: Vec<usize> = (0..n).filter(|&k| se.values[k] > 1e-8).collect();
    let nmo = keep.len();
    if nmo < nalpha {
        return Err(Error::new(
            "Basis set is too small for this many electrons.",
        ));
    }
    let mut x = vec![0.0; n * nmo];
    for i in 0..n {
        for (c, &k) in keep.iter().enumerate() {
            x[i * nmo + c] = se.vectors[i * n + k] / se.values[k].sqrt();
        }
    }
    let xt = transpose(&x, n, nmo);
    let diagonalize = |f: &[f64]| -> (Vec<f64>, Vec<f64>) {
        let fp = matmul(&xt, &matmul(f, &x, n, n, nmo), nmo, n, nmo);
        let e = eigh(&fp, nmo);
        (matmul(&x, &e.vectors, n, nmo, nmo), e.values)
    };
    let fock = |da: &[f64], db: &[f64]| -> (Vec<f64>, Vec<f64>, f64, f64) {
        let dt: Vec<f64> = da.iter().zip(db).map(|(a, b)| a + b).collect();
        let kd: Vec<&[f64]> = if kx == 0.0 {
            vec![]
        } else if unrestricted {
            vec![da, db]
        } else {
            vec![da]
        };
        let (j, k) = build_jk(&eri, &dt, &kd);
        let (vxc, exc) = match (xc, &grid) {
            (Some(f), Some(g)) => {
                let r = xc_matrix(&basis, g, &dt, f);
                (Some(r.v), r.exc)
            }
            _ => (None, 0.0),
        };
        let mut fa = vec![0.0; n * n];
        for q in 0..n * n {
            fa[q] = h[q] + j[q];
            if kx != 0.0 {
                fa[q] -= kx * k[0][q];
            }
            if let Some(v) = &vxc {
                fa[q] += v[q];
            }
        }
        let fb = if unrestricted {
            (0..n * n).map(|q| h[q] + j[q] - kx * k[1][q]).collect()
        } else {
            fa.clone()
        };
        // E = Σ Dt H + ½ Σ Dt J - (kx/2) Σ (Da Ka + Db Kb) + Exc + Enuc
        let mut e = 0.0;
        for q in 0..n * n {
            e += dt[q] * h[q] + 0.5 * dt[q] * j[q];
            if kx != 0.0 {
                let kb = if unrestricted { k[1][q] } else { k[0][q] };
                e -= 0.5 * kx * (da[q] * k[0][q] + db[q] * kb);
            }
        }
        (fa, fb, e + exc + enuc, exc)
    };

    let (mut ca, mut cb, mut eps_a, mut eps_b);
    match guess {
        Some(g) if g.basis == opts.basis && g.nbf == n && g.nmo == nmo => {
            ca = orthonormalize(&g.ca, &s, n, nmo);
            cb = orthonormalize(&g.cb, &s, n, nmo);
            eps_a = vec![0.0; nmo];
            eps_b = vec![0.0; nmo];
        }
        _ => {
            let dt = sad_density(mol, &basis)?;
            let half: Vec<f64> = dt.iter().map(|v| 0.5 * v).collect();
            let (fa, _, _, _) = fock(&half, &half);
            let (c, e) = diagonalize(&fa);
            ca = c.clone();
            cb = c;
            eps_a = e.clone();
            eps_b = e;
        }
    }
    let mut da = density_from(&ca, n, nmo, nalpha);
    let mut db = if unrestricted {
        density_from(&cb, n, nmo, nbeta)
    } else {
        da.clone()
    };

    let mut history = Vec::new();
    let mut diis: Vec<DiisEntry> = Vec::new();
    let mut energy = 0.0;
    let mut prev_energy = 0.0;
    let mut converged = false;
    let mut iterations = 0;
    let mut exc = 0.0;
    let nspin = if unrestricted { 2 } else { 1 };
    for iter in 1..=opts.max_iterations {
        iterations = iter;
        let (fa, fb, e, ex) = fock(&da, &db);
        energy = e;
        exc = ex;
        let fs = if unrestricted { vec![fa, fb] } else { vec![fa] };
        let ds = if unrestricted { [&da, &db] } else { [&da, &da] };
        let mut err = Vec::with_capacity(nspin * nmo * nmo);
        let mut max_err: f64 = 0.0;
        for sp in 0..nspin {
            let fds = matmul(&matmul(&fs[sp], ds[sp], n, n, n), &s, n, n, n);
            let mut comm = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    comm[i * n + j] = fds[i * n + j] - fds[j * n + i];
                }
            }
            let ep = matmul(&xt, &matmul(&comm, &x, n, n, nmo), nmo, n, nmo);
            for v in ep {
                max_err = max_err.max(v.abs());
                err.push(v);
            }
        }
        if opts.record_history {
            history.push(ScfIteration {
                energy,
                error: max_err,
                ca: ca.clone(),
                cb: cb.clone(),
            });
        }
        if iter > 1 && (energy - prev_energy).abs() < opts.conv_energy && max_err < opts.conv_error
        {
            converged = true;
            break;
        }
        prev_energy = energy;
        diis.push(DiisEntry {
            f: fs.clone(),
            e: err,
        });
        if diis.len() > 8 {
            diis.remove(0);
        }
        let fx = if diis.len() >= 2 {
            diis_extrapolate(&diis)
        } else {
            fs
        };
        let (c, e) = diagonalize(&fx[0]);
        ca = c;
        eps_a = e;
        if unrestricted {
            let (c, e) = diagonalize(&fx[1]);
            cb = c;
            eps_b = e;
        } else {
            cb = ca.clone();
            eps_b = eps_a.clone();
        }
        da = density_from(&ca, n, nmo, nalpha);
        db = if unrestricted {
            density_from(&cb, n, nmo, nbeta)
        } else {
            da.clone()
        };
    }

    let dip = dipole_matrices(&basis);
    let mut dipole = [0.0; 3];
    for c in 0..3 {
        let mut acc: f64 = mol.atoms.iter().map(|a| a.z as f64 * a.xyz[c]).sum();
        for q in 0..n * n {
            acc -= (da[q] + db[q]) * dip[c][q];
        }
        dipole[c] = acc;
    }
    Ok(ScfResult {
        converged,
        energy,
        iterations,
        unrestricted,
        method: opts.method,
        basis,
        nbf: n,
        nmo,
        nalpha,
        nbeta,
        ca,
        cb,
        eps_a,
        eps_b,
        da,
        db,
        nuclear_repulsion: enuc,
        exc,
        dipole,
        history,
        grid,
    })
}
