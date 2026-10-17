//! Two-electron repulsion integrals by the McMurchie-Davidson scheme:
//!   (ab|cd) = 2 pi^(5/2) / (p q sqrt(p+q)) Σ_tuv E^ab_tuv Σ_τνφ (-1)^(τ+ν+φ) E^cd_τνφ R_{t+τ,u+ν,v+φ}

use crate::basis::{cartesian_components, component_norms, Basis, Shell};
use crate::hermite::{HermiteE, HermiteR};

pub struct PrimPairData {
    pub a: f64,
    pub b: f64,
    pub p: f64,
    pub pc: [f64; 3],
    pub cc: f64,
    pub e: [HermiteE; 3],
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
                prims.push(PrimPairData { a, b, p, pc, cc, e });
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
    y: Vec<f64>,
}

impl EriWork {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Raw integral block (no component normalization): out[ab * ncd + cd].
pub fn eri_block(ab: &ShellPair, cd: &ShellPair, out: &mut Vec<f64>, wk: &mut EriWork) {
    let ca = cartesian_components(ab.la);
    let cb = cartesian_components(ab.lb);
    let cc = cartesian_components(cd.la);
    let cdc = cartesian_components(cd.lb);
    let nab = ca.len() * cb.len();
    let ncd = cc.len() * cdc.len();
    out.clear();
    out.resize(nab * ncd, 0.0);
    let lab = ab.la + ab.lb;
    let lcd = cd.la + cd.lb;
    let l = lab + lcd;
    let h = lab + 1;
    let cube = h * h * h;
    if wk.w.len() < ncd * cube {
        wk.w.resize(ncd * cube, 0.0);
    }
    for p1 in &ab.prims {
        for p2 in &cd.prims {
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
            let mut icd = 0;
            for c in cc {
                for d in cdc {
                    let base = icd * cube;
                    let (ox, oy, oz) = (
                        p2.e[0].offset(c[0], d[0]),
                        p2.e[1].offset(c[1], d[1]),
                        p2.e[2].offset(c[2], d[2]),
                    );
                    for t in 0..=lab {
                        for u in 0..=lab - t {
                            for v in 0..=lab - t - u {
                                let mut s = 0.0;
                                for tau in 0..=c[0] + d[0] {
                                    let ex = p2.e[0].data[ox + tau];
                                    if ex == 0.0 {
                                        continue;
                                    }
                                    for nu in 0..=c[1] + d[1] {
                                        let exy = ex * p2.e[1].data[oy + nu];
                                        if exy == 0.0 {
                                            continue;
                                        }
                                        for phi in 0..=c[2] + d[2] {
                                            let term = exy
                                                * p2.e[2].data[oz + phi]
                                                * r.at(t + tau, u + nu, v + phi);
                                            if (tau + nu + phi) & 1 == 1 {
                                                s -= term;
                                            } else {
                                                s += term;
                                            }
                                        }
                                    }
                                }
                                w[base + (t * h + u) * h + v] = s;
                            }
                        }
                    }
                    icd += 1;
                }
            }
            let mut iab = 0;
            for a in ca {
                for b in cb {
                    let (ox, oy, oz) = (
                        p1.e[0].offset(a[0], b[0]),
                        p1.e[1].offset(a[1], b[1]),
                        p1.e[2].offset(a[2], b[2]),
                    );
                    let row = iab * ncd;
                    for t in 0..=a[0] + b[0] {
                        let ex = p1.e[0].data[ox + t];
                        if ex == 0.0 {
                            continue;
                        }
                        for u in 0..=a[1] + b[1] {
                            let exy = ex * p1.e[1].data[oy + u];
                            if exy == 0.0 {
                                continue;
                            }
                            for v in 0..=a[2] + b[2] {
                                let e = pref * exy * p1.e[2].data[oz + v];
                                if e == 0.0 {
                                    continue;
                                }
                                let k = (t * h + u) * h + v;
                                for c in 0..ncd {
                                    out[row + c] += e * w[c * cube + k];
                                }
                            }
                        }
                    }
                    iab += 1;
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
    let schwarz = pairs
        .iter()
        .map(|pr| {
            let (na, nb) = pr.ncomp();
            let nab = na * nb;
            eri_block(pr, pr, &mut blk, &mut wk);
            let fa = component_norms(pr.la);
            let fb = component_norms(pr.lb);
            let mut m: f64 = 0.0;
            for x in 0..nab {
                let f = fa[x / nb] * fb[x % nb];
                m = m.max(blk[x * nab + x].abs() * f * f);
            }
            m.sqrt()
        })
        .collect();
    ShellPairList { pairs, schwarz }
}

pub fn compute_eris(basis: &Basis) -> EriTensor {
    let list = shell_pairs(basis, 0);
    let mut eri = EriTensor::new(basis.nbf);
    let mut wk = EriWork::new();
    let mut blk = Vec::new();
    for x in 0..list.pairs.len() {
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
                            let k =
                                eri.index(ab.off_a + a, ab.off_b + b, cd.off_a + c, cd.off_b + d);
                            eri.data[k] = blk[idx] * fab * fc[c] * fd[d];
                            idx += 1;
                        }
                    }
                }
            }
        }
    }
    eri
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
    for x in 0..list.pairs.len() {
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
    for p1 in &ab.prims {
        for p2 in &cd.prims {
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
            let mut icd = 0;
            for c in cc {
                for d in cdc {
                    let base = icd * cube;
                    let (ox, oy, oz) = (
                        p2.e[0].offset(c[0], d[0]),
                        p2.e[1].offset(c[1], d[1]),
                        p2.e[2].offset(c[2], d[2]),
                    );
                    for t in 0..=lab + 1 {
                        for u in 0..=lab + 1 - t {
                            for v in 0..=lab + 1 - t - u {
                                let mut s = 0.0;
                                for tau in 0..=c[0] + d[0] {
                                    let ex = p2.e[0].data[ox + tau];
                                    if ex == 0.0 {
                                        continue;
                                    }
                                    for nu in 0..=c[1] + d[1] {
                                        let exy = ex * p2.e[1].data[oy + nu];
                                        if exy == 0.0 {
                                            continue;
                                        }
                                        for phi in 0..=c[2] + d[2] {
                                            let term = exy
                                                * p2.e[2].data[oz + phi]
                                                * r.at(t + tau, u + nu, v + phi);
                                            if (tau + nu + phi) & 1 == 1 {
                                                s -= term;
                                            } else {
                                                s += term;
                                            }
                                        }
                                    }
                                }
                                w[base + (t * h + u) * h + v] = s;
                            }
                        }
                    }
                    icd += 1;
                }
            }
            let mut iab = 0;
            for a in ca {
                for b in cb {
                    let yb = iab * cube;
                    y[yb..yb + cube].fill(0.0);
                    let mut any = false;
                    for c in 0..ncd {
                        let gv = gamma[iab * ncd + c];
                        if gv == 0.0 {
                            continue;
                        }
                        any = true;
                        let wb = c * cube;
                        for k in 0..cube {
                            y[yb + k] += gv * w[wb + k];
                        }
                    }
                    if any {
                        let (ea, eb) = (p1.a, p1.b);
                        let (exs, eys, ezs) = (&p1.e[0], &p1.e[1], &p1.e[2]);
                        for t in 0..=a[0] + b[0] + 1 {
                            let ex = exs.get(a[0] as isize, b[0] as isize, t);
                            let dax = exs.d_a(a[0], b[0], t, ea);
                            let dbx = exs.d_b(a[0], b[0], t, eb);
                            for u in 0..=a[1] + b[1] + 1 {
                                let ey = eys.get(a[1] as isize, b[1] as isize, u);
                                let day = eys.d_a(a[1], b[1], u, ea);
                                let dby = eys.d_b(a[1], b[1], u, eb);
                                for v in 0..=a[2] + b[2] + 1 {
                                    let yv = pref * y[yb + (t * h + u) * h + v];
                                    if yv == 0.0 {
                                        continue;
                                    }
                                    let ez = ezs.get(a[2] as isize, b[2] as isize, v);
                                    let daz = ezs.d_a(a[2], b[2], v, ea);
                                    let dbz = ezs.d_b(a[2], b[2], v, eb);
                                    g[0] += dax * ey * ez * yv;
                                    g[1] += ex * day * ez * yv;
                                    g[2] += ex * ey * daz * yv;
                                    g[3] += dbx * ey * ez * yv;
                                    g[4] += ex * dby * ez * yv;
                                    g[5] += ex * ey * dbz * yv;
                                }
                            }
                        }
                    }
                    iab += 1;
                }
            }

            // Ket derivative: G^cd_tuv = Σ_ab Γ_abcd E^ab_tuv (t+u+v <= lab), stored in y.
            y[..ncd * cube].fill(0.0);
            iab = 0;
            for a in ca {
                for b in cb {
                    let (ox, oy, oz) = (
                        p1.e[0].offset(a[0], b[0]),
                        p1.e[1].offset(a[1], b[1]),
                        p1.e[2].offset(a[2], b[2]),
                    );
                    for t in 0..=a[0] + b[0] {
                        let ex = p1.e[0].data[ox + t];
                        if ex == 0.0 {
                            continue;
                        }
                        for u in 0..=a[1] + b[1] {
                            let exy = ex * p1.e[1].data[oy + u];
                            if exy == 0.0 {
                                continue;
                            }
                            for v in 0..=a[2] + b[2] {
                                let e = exy * p1.e[2].data[oz + v];
                                if e == 0.0 {
                                    continue;
                                }
                                let k = (t * h + u) * h + v;
                                for c in 0..ncd {
                                    y[c * cube + k] += gamma[iab * ncd + c] * e;
                                }
                            }
                        }
                    }
                    iab += 1;
                }
            }
            icd = 0;
            for c in cc {
                for d in cdc {
                    let yb = icd * cube;
                    let ec = p2.a;
                    let (exs, eys, ezs) = (&p2.e[0], &p2.e[1], &p2.e[2]);
                    for tau in 0..=c[0] + d[0] + 1 {
                        let ex = exs.get(c[0] as isize, d[0] as isize, tau);
                        let dcx = exs.d_a(c[0], d[0], tau, ec);
                        for nu in 0..=c[1] + d[1] + 1 {
                            let ey = eys.get(c[1] as isize, d[1] as isize, nu);
                            let dcy = eys.d_a(c[1], d[1], nu, ec);
                            for phi in 0..=c[2] + d[2] + 1 {
                                let ez = ezs.get(c[2] as isize, d[2] as isize, phi);
                                let dcz = ezs.d_a(c[2], d[2], phi, ec);
                                let w0 = dcx * ey * ez;
                                let w1 = ex * dcy * ez;
                                let w2 = ex * ey * dcz;
                                if w0 == 0.0 && w1 == 0.0 && w2 == 0.0 {
                                    continue;
                                }
                                let mut s = 0.0;
                                for t in 0..=lab {
                                    for u in 0..=lab - t {
                                        let gb = yb + (t * h + u) * h;
                                        for v in 0..=lab - t - u {
                                            s += y[gb + v] * r.at(t + tau, u + nu, v + phi);
                                        }
                                    }
                                }
                                if (tau + nu + phi) & 1 == 1 {
                                    s = -s;
                                }
                                s *= pref;
                                g[6] += w0 * s;
                                g[7] += w1 * s;
                                g[8] += w2 * s;
                            }
                        }
                    }
                    icd += 1;
                }
            }
        }
    }
    g
}
