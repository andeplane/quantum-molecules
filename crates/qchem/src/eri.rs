//! Two-electron repulsion integrals by the McMurchie-Davidson scheme:
//!   (ab|cd) = 2 pi^(5/2) / (p q sqrt(p+q)) Σ_tuv E^ab_tuv Σ_τνφ (-1)^(τ+ν+φ) E^cd_τνφ R_{t+τ,u+ν,v+φ}

use crate::basis::{cartesian_components, component_norms, Basis, Shell};
use crate::hermite::{HermiteE, HermiteR};

#[derive(Clone)]
pub struct PrimPairData {
    /// Schwarz bound sqrt(max |(ab|ab)|) of this primitive pair alone (infinite until
    /// `shell_pairs` fills it in), used to skip negligible primitive quartets.
    pub q: f64,
    pub a: f64,
    pub b: f64,
    pub p: f64,
    pub pc: [f64; 3],
    pub cc: f64,
    pub e: [HermiteE; 3],
    /// Non-zero products E_t E_u E_v for every component pair (a, b), flattened; component pair
    /// `k` owns `terms[ranges[k].0..ranges[k].1]`. Saves the zero tests and nested loops in the
    /// integral kernel.
    pub terms: Vec<HermiteTerm>,
    pub ranges: Vec<(u32, u32)>,
    /// Derivatives of the products for gradients (only for pairs built with `extra >= 1`):
    /// component pair `k` owns `dterms[dranges[k].0..dranges[k].1]`.
    pub dterms: Vec<DerivTerm>,
    pub dranges: Vec<(u32, u32)>,
}

/// d/dA and d/dB (x, y, z) of E^{ab}_t E^{ab}_u E^{ab}_v.
#[derive(Clone, Copy, Debug)]
pub struct DerivTerm {
    pub t: u8,
    pub u: u8,
    pub v: u8,
    pub da: [f64; 3],
    pub db: [f64; 3],
}

fn deriv_terms(
    e: &[HermiteE; 3],
    la: usize,
    lb: usize,
    ea: f64,
    eb: f64,
) -> (Vec<DerivTerm>, Vec<(u32, u32)>) {
    let mut terms = Vec::new();
    let mut ranges = Vec::new();
    for a in cartesian_components(la) {
        for b in cartesian_components(lb) {
            let start = terms.len() as u32;
            for t in 0..=a[0] + b[0] + 1 {
                let ex = e[0].get(a[0] as isize, b[0] as isize, t);
                let (dax, dbx) = (e[0].d_a(a[0], b[0], t, ea), e[0].d_b(a[0], b[0], t, eb));
                for u in 0..=a[1] + b[1] + 1 {
                    let ey = e[1].get(a[1] as isize, b[1] as isize, u);
                    let (day, dby) = (e[1].d_a(a[1], b[1], u, ea), e[1].d_b(a[1], b[1], u, eb));
                    for v in 0..=a[2] + b[2] + 1 {
                        let ez = e[2].get(a[2] as isize, b[2] as isize, v);
                        let (daz, dbz) = (e[2].d_a(a[2], b[2], v, ea), e[2].d_b(a[2], b[2], v, eb));
                        let da = [dax * ey * ez, ex * day * ez, ex * ey * daz];
                        let db = [dbx * ey * ez, ex * dby * ez, ex * ey * dbz];
                        if da.iter().chain(&db).all(|&x| x == 0.0) {
                            continue;
                        }
                        terms.push(DerivTerm {
                            t: t as u8,
                            u: u as u8,
                            v: v as u8,
                            da,
                            db,
                        });
                    }
                }
            }
            ranges.push((start, terms.len() as u32));
        }
    }
    (terms, ranges)
}

/// One product of Hermite expansion coefficients, E^{ab}_t E^{ab}_u E^{ab}_v.
#[derive(Clone, Copy, Debug)]
pub struct HermiteTerm {
    pub t: u8,
    pub u: u8,
    pub v: u8,
    pub e: f64,
    /// `e` times (-1)^(t+u+v), the form used on the ket side.
    pub signed: f64,
}

fn hermite_terms(e: &[HermiteE; 3], la: usize, lb: usize) -> (Vec<HermiteTerm>, Vec<(u32, u32)>) {
    let mut terms = Vec::new();
    let mut ranges = Vec::new();
    for a in cartesian_components(la) {
        for b in cartesian_components(lb) {
            let start = terms.len() as u32;
            let (ox, oy, oz) = (
                e[0].offset(a[0], b[0]),
                e[1].offset(a[1], b[1]),
                e[2].offset(a[2], b[2]),
            );
            for t in 0..=a[0] + b[0] {
                let ex = e[0].data[ox + t];
                for u in 0..=a[1] + b[1] {
                    let exy = ex * e[1].data[oy + u];
                    for v in 0..=a[2] + b[2] {
                        let val = exy * e[2].data[oz + v];
                        if val == 0.0 {
                            continue;
                        }
                        let signed = if (t + u + v) & 1 == 1 { -val } else { val };
                        terms.push(HermiteTerm {
                            t: t as u8,
                            u: u as u8,
                            v: v as u8,
                            e: val,
                            signed,
                        });
                    }
                }
            }
            ranges.push((start, terms.len() as u32));
        }
    }
    (terms, ranges)
}

/// Precomputed data for a shell pair (i >= j).
pub struct ShellPair {
    pub i: usize,
    pub j: usize,
    pub la: usize,
    pub lb: usize,
    pub atom_a: usize,
    pub atom_b: usize,
    pub off_a: usize,
    pub off_b: usize,
    pub prims: Vec<PrimPairData>,
}

impl ShellPair {
    /// `extra` raises the Hermite tables by that many powers on both centres (1 for gradients).
    pub fn new(sa: &Shell, sb: &Shell, i: usize, j: usize, extra: usize) -> Self {
        let mut prims = Vec::new();
        for (k, &a) in sa.exps.iter().enumerate() {
            for (m, &b) in sb.exps.iter().enumerate() {
                let p = a + b;
                let pc = [0, 1, 2].map(|d| (a * sa.center[d] + b * sb.center[d]) / p);
                let e = [0, 1, 2].map(|d| {
                    HermiteE::new(
                        sa.l + extra,
                        sb.l + extra,
                        a,
                        b,
                        sa.center[d] - sb.center[d],
                    )
                });
                let cc = sa.coefs[k] * sb.coefs[m];
                if (cc * e[0].data[0] * e[1].data[0] * e[2].data[0]).abs() < 1e-15 {
                    continue;
                }
                let (terms, ranges) = hermite_terms(&e, sa.l, sb.l);
                let (dterms, dranges) = if extra >= 1 {
                    deriv_terms(&e, sa.l, sb.l, a, b)
                } else {
                    (Vec::new(), Vec::new())
                };
                prims.push(PrimPairData {
                    q: f64::INFINITY,
                    a,
                    b,
                    p,
                    pc,
                    cc,
                    e,
                    terms,
                    ranges,
                    dterms,
                    dranges,
                });
            }
        }
        ShellPair {
            i,
            j,
            la: sa.l,
            lb: sb.l,
            atom_a: sa.atom,
            atom_b: sb.atom,
            off_a: sa.offset,
            off_b: sb.offset,
            prims,
        }
    }

    pub fn ncomp(&self) -> (usize, usize) {
        (
            cartesian_components(self.la).len(),
            cartesian_components(self.lb).len(),
        )
    }
}

const TWO_PI_2_5: f64 = 2.0 * 17.493_418_327_624_862; // 2 * pi^(5/2)

/// Scratch space reused across quartets.
#[derive(Default)]
pub struct EriWork {
    r: HermiteR,
    w: Vec<f64>,
    h_off: Vec<usize>,
    h_pos: Vec<usize>,
    /// Scratch for the gradient kernel.
    y: Vec<f64>,
}

impl EriWork {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Raw integral block (no component normalization): out[ab * ncd + cd].
pub fn eri_block(ab: &ShellPair, cd: &ShellPair, out: &mut Vec<f64>, wk: &mut EriWork) {
    let nab = cartesian_components(ab.la).len() * cartesian_components(ab.lb).len();
    let ncd = cartesian_components(cd.la).len() * cartesian_components(cd.lb).len();
    out.clear();
    out.resize(nab * ncd, 0.0);
    let lab = ab.la + ab.lb;
    let lcd = cd.la + cd.lb;
    let l = lab + lcd;
    let l1 = l + 1;
    // Bra Hermite indices (t, u, v), t + u + v <= lab: their offsets into the R cube (the ket
    // indices add to these) and their positions in the rows of `w`.
    let h = lab + 1;
    wk.h_off.clear();
    wk.h_pos.clear();
    wk.h_pos.resize(h * h * h, 0);
    for t in 0..=lab {
        for u in 0..=lab - t {
            for v in 0..=lab - t - u {
                wk.h_pos[(t * h + u) * h + v] = wk.h_off.len();
                wk.h_off.push((t * l1 + u) * l1 + v);
            }
        }
    }
    let nh = wk.h_off.len();
    if wk.w.len() < ncd * nh {
        wk.w.resize(ncd * nh, 0.0);
    }
    for p1 in &ab.prims {
        for p2 in &cd.prims {
            if p1.q * p2.q < PRIM_CUTOFF {
                continue;
            }
            let (p, q) = (p1.p, p2.p);
            let alpha = p * q / (p + q);
            let pref = TWO_PI_2_5 / (p * q * (p + q).sqrt()) * p1.cc * p2.cc;
            wk.r.compute(
                l,
                alpha,
                p1.pc[0] - p2.pc[0],
                p1.pc[1] - p2.pc[1],
                p1.pc[2] - p2.pc[2],
            );
            let r = &wk.r.out[..l1 * l1 * l1];
            let w = &mut wk.w[..ncd * nh];
            let (h_off, h_pos) = (&wk.h_off, &wk.h_pos);
            // Half transform over the ket: w[cd][tuv] = Σ_τνφ (-1)^(τ+ν+φ) E^cd_τνφ R_{t+τ,u+ν,v+φ}
            for (icd, &(s0, s1)) in p2.ranges.iter().enumerate() {
                let terms = &p2.terms[s0 as usize..s1 as usize];
                let row = &mut w[icd * nh..(icd + 1) * nh];
                for (slot, &off) in row.iter_mut().zip(h_off) {
                    let mut acc = 0.0;
                    for term in terms {
                        let k =
                            off + (term.t as usize * l1 + term.u as usize) * l1 + term.v as usize;
                        acc += term.signed * r[k];
                    }
                    *slot = acc;
                }
            }
            for (iab, &(s0, s1)) in p1.ranges.iter().enumerate() {
                let dst = &mut out[iab * ncd..(iab + 1) * ncd];
                for term in &p1.terms[s0 as usize..s1 as usize] {
                    let e = pref * term.e;
                    let k = h_pos[(term.t as usize * h + term.u as usize) * h + term.v as usize];
                    for (c, o) in dst.iter_mut().enumerate() {
                        *o += e * w[c * nh + k];
                    }
                }
            }
        }
    }
}

/// Packed 8-fold symmetric ERI storage.
pub struct EriTensor {
    pub n: usize,
    pub data: Vec<f64>,
}

#[inline(always)]
pub fn pair_index(i: usize, j: usize) -> usize {
    if i >= j {
        i * (i + 1) / 2 + j
    } else {
        j * (j + 1) / 2 + i
    }
}

impl EriTensor {
    pub fn new(n: usize) -> Self {
        let np = n * (n + 1) / 2;
        EriTensor {
            n,
            data: vec![0.0; np * (np + 1) / 2],
        }
    }

    /// Assemble the tensor from the stripes of `compute_eri_part`.
    pub fn from_parts<'a>(
        n: usize,
        parts: impl IntoIterator<Item = (&'a [u32], &'a [f64])>,
    ) -> Self {
        let mut eri = EriTensor::new(n);
        for (index, value) in parts {
            eri.add_part(index, value);
        }
        eri
    }

    /// Write one stripe from `compute_eri_part` into the tensor.
    pub fn add_part(&mut self, index: &[u32], value: &[f64]) {
        for (&k, &v) in index.iter().zip(value) {
            self.data[k as usize] = v;
        }
    }

    #[inline(always)]
    pub fn index(&self, i: usize, j: usize, k: usize, l: usize) -> usize {
        pair_index(pair_index(i, j), pair_index(k, l))
    }

    pub fn get(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        self.data[self.index(i, j, k, l)]
    }
}

pub struct ShellPairList {
    pub pairs: Vec<ShellPair>,
    /// Schwarz bound sqrt(max |(ab|ab)|) per pair.
    pub schwarz: Vec<f64>,
}

pub const SCHWARZ_CUTOFF: f64 = 1e-12;
/// Primitive quartets below this Schwarz bound are skipped; far below the accuracy that matters.
const PRIM_CUTOFF: f64 = 1e-15;

pub fn shell_pairs(basis: &Basis, extra: usize) -> ShellPairList {
    let mut pairs = Vec::new();
    for i in 0..basis.shells.len() {
        for j in 0..=i {
            pairs.push(ShellPair::new(
                &basis.shells[i],
                &basis.shells[j],
                i,
                j,
                extra,
            ));
        }
    }
    let mut wk = EriWork::new();
    let mut blk = Vec::new();
    let schwarz: Vec<f64> = pairs
        .iter()
        .map(|pr| schwarz_bound(pr, &mut blk, &mut wk))
        .collect();
    // Primitive bounds come after the shell bounds, which must not be screened themselves.
    // For gradient pairs (extra = 1) the bound is on the integral, not its derivative, which
    // only adds a factor of order the exponent; still far below what matters.
    {
        for pr in pairs.iter_mut() {
            for k in 0..pr.prims.len() {
                let single = ShellPair {
                    prims: vec![pr.prims[k].clone()],
                    ..*pr
                };
                pr.prims[k].q = schwarz_bound(&single, &mut blk, &mut wk);
            }
        }
    }
    ShellPairList { pairs, schwarz }
}

/// sqrt(max |(ab|ab)|) over the components of a shell pair, with component norms applied.
fn schwarz_bound(pr: &ShellPair, blk: &mut Vec<f64>, wk: &mut EriWork) -> f64 {
    let (na, nb) = pr.ncomp();
    let nab = na * nb;
    eri_block(pr, pr, blk, wk);
    let fa = component_norms(pr.la);
    let fb = component_norms(pr.lb);
    let mut m: f64 = 0.0;
    for x in 0..nab {
        let f = fa[x / nb] * fb[x % nb];
        m = m.max(blk[x * nab + x].abs() * f * f);
    }
    m.sqrt()
}

pub fn compute_eris(basis: &Basis) -> EriTensor {
    let mut eri = EriTensor::new(basis.nbf);
    eri_stripe(basis, 0, 1, |k, v| eri.data[k] = v);
    eri
}

/// The integrals of every `stride`-th shell-pair row starting at `part`, as packed indices and
/// values. The stripes for `part = 0..stride` together hold exactly what `compute_eris` computes,
/// so separate workers can share the work and the results be merged with `EriTensor::from_parts`.
pub fn compute_eri_part(basis: &Basis, part: usize, stride: usize) -> (Vec<u32>, Vec<f64>) {
    let mut index = Vec::new();
    let mut value = Vec::new();
    eri_stripe(basis, part, stride, |k, v| {
        if v != 0.0 {
            index.push(k as u32);
            value.push(v);
        }
    });
    (index, value)
}

fn eri_stripe(basis: &Basis, part: usize, stride: usize, mut store: impl FnMut(usize, f64)) {
    let list = shell_pairs(basis, 0);
    let mut wk = EriWork::new();
    let mut blk = Vec::new();
    for x in (part..list.pairs.len()).step_by(stride.max(1)) {
        for y in 0..=x {
            if list.schwarz[x] * list.schwarz[y] < SCHWARZ_CUTOFF {
                continue;
            }
            let ab = &list.pairs[x];
            let cd = &list.pairs[y];
            eri_block(ab, cd, &mut blk, &mut wk);
            let (na, nb) = ab.ncomp();
            let (nc, nd) = cd.ncomp();
            let (fa, fb, fc, fd) = (
                component_norms(ab.la),
                component_norms(ab.lb),
                component_norms(cd.la),
                component_norms(cd.lb),
            );
            let mut idx = 0;
            for a in 0..na {
                for b in 0..nb {
                    let fab = fa[a] * fb[b];
                    for c in 0..nc {
                        for d in 0..nd {
                            let k = pair_index(
                                pair_index(ab.off_a + a, ab.off_b + b),
                                pair_index(cd.off_a + c, cd.off_b + d),
                            );
                            store(k, blk[idx] * fab * fc[c] * fd[d]);
                            idx += 1;
                        }
                    }
                }
            }
        }
    }
}

/// Coulomb matrix J[D] and exchange matrices K[D_m] (K_ik = Σ (ij|kl) D_jl) from packed ERIs.
pub fn build_jk(eri: &EriTensor, jd: &[f64], kds: &[&[f64]]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = eri.n;
    let mut j_mat = vec![0.0; n * n];
    let mut ks: Vec<Vec<f64>> = kds.iter().map(|_| vec![0.0; n * n]).collect();
    let data = &eri.data;
    let mut idx = 0;
    let mut ij = 0;
    for i in 0..n {
        for j in 0..=i {
            let mut kl = 0;
            'outer: for k in 0..=i {
                for l in 0..=k {
                    if kl > ij {
                        break 'outer;
                    }
                    let mut v = data[idx];
                    idx += 1;
                    if v != 0.0 {
                        if i == j {
                            v *= 0.5;
                        }
                        if k == l {
                            v *= 0.5;
                        }
                        if ij == kl {
                            v *= 0.5;
                        }
                        let jv = 2.0 * v;
                        let dkl = jd[k * n + l] * jv;
                        let dij = jd[i * n + j] * jv;
                        j_mat[i * n + j] += dkl;
                        j_mat[j * n + i] += dkl;
                        j_mat[k * n + l] += dij;
                        j_mat[l * n + k] += dij;
                        for (km, d) in ks.iter_mut().zip(kds.iter()) {
                            km[i * n + k] += d[j * n + l] * v;
                            km[j * n + k] += d[i * n + l] * v;
                            km[i * n + l] += d[j * n + k] * v;
                            km[j * n + l] += d[i * n + k] * v;
                            km[k * n + i] += d[l * n + j] * v;
                            km[l * n + i] += d[k * n + j] * v;
                            km[k * n + j] += d[l * n + i] * v;
                            km[l * n + j] += d[k * n + i] * v;
                        }
                    }
                    kl += 1;
                }
            }
            ij += 1;
        }
    }
    (j_mat, ks)
}

/// Two-electron gradient ½ Σ Γ_μνλσ d(μν|λσ)/dX with
/// Γ = Dt_μν Dt_λσ - (kx/2)(Da_μλ Da_νσ + Da_μσ Da_νλ + Db_μλ Db_νσ + Db_μσ Db_νλ).
pub fn eri_gradient(
    basis: &Basis,
    natom: usize,
    dt: &[f64],
    da: &[f64],
    db: &[f64],
    kx: f64,
) -> Vec<f64> {
    eri_gradient_part(basis, natom, dt, da, db, kx, 0, 1)
}

/// The share of `eri_gradient` from every `stride`-th shell-pair row starting at `part`; the
/// shares for `part = 0..stride` add up to the full two-electron gradient.
#[allow(clippy::too_many_arguments)]
pub fn eri_gradient_part(
    basis: &Basis,
    natom: usize,
    dt: &[f64],
    da: &[f64],
    db: &[f64],
    kx: f64,
    part: usize,
    stride: usize,
) -> Vec<f64> {
    let n = basis.nbf;
    let mut grad = vec![0.0; natom * 3];
    let list = shell_pairs(basis, 1);
    let ns = basis.shells.len();
    let mut dm = vec![0.0; ns * ns];
    for (i, sa) in basis.shells.iter().enumerate() {
        for (j, sb) in basis.shells.iter().enumerate() {
            let mut m: f64 = 0.0;
            for a in 0..cartesian_components(sa.l).len() {
                for b in 0..cartesian_components(sb.l).len() {
                    let k = (sa.offset + a) * n + sb.offset + b;
                    m = m.max(dt[k].abs()).max(da[k].abs()).max(db[k].abs());
                }
            }
            dm[i * ns + j] = m;
        }
    }
    let mut wk = EriWork::new();
    let mut gamma = Vec::new();
    for x in (part..list.pairs.len()).step_by(stride.max(1)) {
        for y in 0..=x {
            let sch = list.schwarz[x] * list.schwarz[y];
            if sch < SCHWARZ_CUTOFF {
                continue;
            }
            let ab = &list.pairs[x];
            let cd = &list.pairs[y];
            let (p, q, r, s) = (ab.i, ab.j, cd.i, cd.j);
            let dbound = (dm[p * ns + q] * dm[r * ns + s])
                .max(dm[p * ns + r] * dm[q * ns + s])
                .max(dm[p * ns + s] * dm[q * ns + r]);
            if sch * dbound < 1e-11 {
                continue;
            }
            let mut deg = 8.0;
            if p == q {
                deg *= 0.5;
            }
            if r == s {
                deg *= 0.5;
            }
            if x == y {
                deg *= 0.5;
            }
            let (na, nb) = ab.ncomp();
            let (nc, nd) = cd.ncomp();
            let (fa, fb, fc, fd) = (
                component_norms(ab.la),
                component_norms(ab.lb),
                component_norms(cd.la),
                component_norms(cd.lb),
            );
            gamma.clear();
            for a in 0..na {
                let mu = ab.off_a + a;
                for b in 0..nb {
                    let nu = ab.off_b + b;
                    for c in 0..nc {
                        let la = cd.off_a + c;
                        for d in 0..nd {
                            let si = cd.off_b + d;
                            let mut g = dt[mu * n + nu] * dt[la * n + si];
                            if kx != 0.0 {
                                g -= 0.5
                                    * kx
                                    * (da[mu * n + la] * da[nu * n + si]
                                        + da[mu * n + si] * da[nu * n + la]
                                        + db[mu * n + la] * db[nu * n + si]
                                        + db[mu * n + si] * db[nu * n + la]);
                            }
                            gamma.push(0.5 * deg * g * fa[a] * fb[b] * fc[c] * fd[d]);
                        }
                    }
                }
            }
            let g = quartet_gradient(ab, cd, &gamma, &mut wk);
            let atoms = [ab.atom_a, ab.atom_b, cd.atom_a];
            for k in 0..3 {
                for (c, &at) in atoms.iter().enumerate() {
                    grad[3 * at + k] += g[3 * c + k];
                }
                grad[3 * cd.atom_b + k] -= g[k] + g[3 + k] + g[6 + k];
            }
        }
    }
    grad
}

/// Σ_abcd Γ_abcd d(ab|cd)/dX for X = A, B, C. Returns [dA(3), dB(3), dC(3)].
fn quartet_gradient(ab: &ShellPair, cd: &ShellPair, gamma: &[f64], wk: &mut EriWork) -> [f64; 9] {
    let mut g = [0.0; 9];
    let ca = cartesian_components(ab.la);
    let cb = cartesian_components(ab.lb);
    let cc = cartesian_components(cd.la);
    let cdc = cartesian_components(cd.lb);
    let ncd = cc.len() * cdc.len();
    let nab = ca.len() * cb.len();
    let lab = ab.la + ab.lb;
    let l = lab + cd.la + cd.lb + 1;
    let h = lab + 2;
    let cube = h * h * h;
    let need = ncd.max(nab) * cube;
    if wk.w.len() < need {
        wk.w.resize(need, 0.0);
    }
    if wk.y.len() < need {
        wk.y.resize(need, 0.0);
    }
    // Bra Hermite indices t + u + v <= lab: their slots in `y` (h_pos) and R offsets (h_off).
    let l1 = l + 1;
    wk.h_off.clear();
    wk.h_pos.clear();
    for t in 0..=lab {
        for u in 0..=lab - t {
            for v in 0..=lab - t - u {
                wk.h_pos.push((t * h + u) * h + v);
                wk.h_off.push((t * l1 + u) * l1 + v);
            }
        }
    }
    for p1 in &ab.prims {
        for p2 in &cd.prims {
            if p1.q * p2.q < PRIM_CUTOFF {
                continue;
            }
            let (p, q) = (p1.p, p2.p);
            let alpha = p * q / (p + q);
            let pref = TWO_PI_2_5 / (p * q * (p + q).sqrt()) * p1.cc * p2.cc;
            wk.r.compute(
                l,
                alpha,
                p1.pc[0] - p2.pc[0],
                p1.pc[1] - p2.pc[1],
                p1.pc[2] - p2.pc[2],
            );
            let r = &wk.r;
            let w = &mut wk.w;
            let y = &mut wk.y;

            // Bra derivatives: W^cd_tuv for t+u+v <= lab+1.
            let rr = &r.out[..l1 * l1 * l1];
            for (icd, &(s0, s1)) in p2.ranges.iter().enumerate() {
                let base = icd * cube;
                let terms = &p2.terms[s0 as usize..s1 as usize];
                for t in 0..=lab + 1 {
                    for u in 0..=lab + 1 - t {
                        for v in 0..=lab + 1 - t - u {
                            let off = (t * l1 + u) * l1 + v;
                            let mut acc = 0.0;
                            for term in terms {
                                let k = off
                                    + (term.t as usize * l1 + term.u as usize) * l1
                                    + term.v as usize;
                                acc += term.signed * rr[k];
                            }
                            w[base + (t * h + u) * h + v] = acc;
                        }
                    }
                }
            }
            for iab in 0..nab {
                let grow = &gamma[iab * ncd..(iab + 1) * ncd];
                if grow.iter().all(|&x| x == 0.0) {
                    continue;
                }
                let (s0, s1) = p1.dranges[iab];
                for term in &p1.dterms[s0 as usize..s1 as usize] {
                    let k = (term.t as usize * h + term.u as usize) * h + term.v as usize;
                    let mut yv = 0.0;
                    for (c, &gv) in grow.iter().enumerate() {
                        yv += gv * w[c * cube + k];
                    }
                    yv *= pref;
                    for x in 0..3 {
                        g[x] += term.da[x] * yv;
                        g[3 + x] += term.db[x] * yv;
                    }
                }
            }

            // Ket derivative: G^cd_tuv = Σ_ab Γ_abcd E^ab_tuv (t+u+v <= lab), stored in y.
            y[..ncd * cube].fill(0.0);
            for (iab, &(s0, s1)) in p1.ranges.iter().enumerate() {
                for term in &p1.terms[s0 as usize..s1 as usize] {
                    let k = (term.t as usize * h + term.u as usize) * h + term.v as usize;
                    for c in 0..ncd {
                        y[c * cube + k] += gamma[iab * ncd + c] * term.e;
                    }
                }
            }
            for (icd, &(s0, s1)) in p2.dranges.iter().enumerate() {
                let yb = icd * cube;
                for term in &p2.dterms[s0 as usize..s1 as usize] {
                    let shift = (term.t as usize * l1 + term.u as usize) * l1 + term.v as usize;
                    let mut acc = 0.0;
                    for (&yo, &ro) in wk.h_pos.iter().zip(&wk.h_off) {
                        acc += y[yb + yo] * rr[ro + shift];
                    }
                    if (term.t + term.u + term.v) & 1 == 1 {
                        acc = -acc;
                    }
                    acc *= pref;
                    for x in 0..3 {
                        g[6 + x] += term.da[x] * acc;
                    }
                }
            }
        }
    }
    g
}
