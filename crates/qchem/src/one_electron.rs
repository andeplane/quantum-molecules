//! One-electron integrals (overlap, kinetic, nuclear attraction, dipole) and their nuclear
//! derivatives, all by the McMurchie-Davidson scheme.

use crate::basis::{cartesian_components, component_norms, Basis, Shell};
use crate::hermite::{HermiteE, HermiteR};
use crate::molecule::Atom;
use std::f64::consts::PI;

struct PrimPair {
    a: f64,
    b: f64,
    p: f64,
    pc: [f64; 3],
    cc: f64,
    e: [HermiteE; 3],
}

fn prim_pairs(sa: &Shell, sb: &Shell, iextra: usize, jextra: usize) -> Vec<PrimPair> {
    let mut out = Vec::with_capacity(sa.exps.len() * sb.exps.len());
    for (i, &a) in sa.exps.iter().enumerate() {
        for (j, &b) in sb.exps.iter().enumerate() {
            let p = a + b;
            let pc = [0, 1, 2].map(|k| (a * sa.center[k] + b * sb.center[k]) / p);
            let e = [0, 1, 2].map(|k| {
                HermiteE::new(
                    sa.l + iextra,
                    sb.l + jextra,
                    a,
                    b,
                    sa.center[k] - sb.center[k],
                )
            });
            out.push(PrimPair {
                a,
                b,
                p,
                pc,
                cc: sa.coefs[i] * sb.coefs[j],
                e,
            });
        }
    }
    out
}

/// 1D overlap without the sqrt(pi/p) factor.
#[inline]
fn s1(e: &HermiteE, i: isize, j: isize) -> f64 {
    e.get(i, j, 0)
}

/// 1D kinetic piece: -2b^2 S(i,j+2) + b(2j+1) S(i,j) - j(j-1)/2 S(i,j-2)
#[inline]
fn t1(e: &HermiteE, b: f64, i: isize, j: isize) -> f64 {
    if i < 0 {
        return 0.0;
    }
    let mut v = -2.0 * b * b * s1(e, i, j + 2) + b * (2 * j + 1) as f64 * s1(e, i, j);
    if j >= 2 {
        v -= 0.5 * (j * (j - 1)) as f64 * s1(e, i, j - 2);
    }
    v
}

fn scatter(m: &mut [f64], blk: &[f64], sa: &Shell, sb: &Shell, n: usize) {
    let fa = component_norms(sa.l);
    let fb = component_norms(sb.l);
    let nb = fb.len();
    for ia in 0..fa.len() {
        for ib in 0..nb {
            let v = blk[ia * nb + ib] * fa[ia] * fb[ib];
            m[(sa.offset + ia) * n + sb.offset + ib] = v;
            m[(sb.offset + ib) * n + sa.offset + ia] = v;
        }
    }
}

pub struct OneElectron {
    pub s: Vec<f64>,
    pub t: Vec<f64>,
    pub v: Vec<f64>,
}

/// Overlap, kinetic and nuclear attraction matrices (n x n, row-major).
pub fn one_electron_matrices(basis: &Basis, atoms: &[Atom]) -> OneElectron {
    let n = basis.nbf;
    let mut s = vec![0.0; n * n];
    let mut t = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut r = HermiteR::new();
    for si in 0..basis.shells.len() {
        for sj in 0..=si {
            let sa = &basis.shells[si];
            let sb = &basis.shells[sj];
            let ca = cartesian_components(sa.l);
            let cb = cartesian_components(sb.l);
            let nb = cb.len();
            let mut sblk = vec![0.0; ca.len() * nb];
            let mut tblk = vec![0.0; ca.len() * nb];
            let mut vblk = vec![0.0; ca.len() * nb];
            let l = sa.l + sb.l;
            for pp in prim_pairs(sa, sb, 0, 2) {
                let sq = (PI / pp.p).powf(1.5) * pp.cc;
                for (ia, a) in ca.iter().enumerate() {
                    for (ib, b) in cb.iter().enumerate() {
                        let sv: [f64; 3] =
                            [0, 1, 2].map(|k| s1(&pp.e[k], a[k] as isize, b[k] as isize));
                        let tv: [f64; 3] =
                            [0, 1, 2].map(|k| t1(&pp.e[k], pp.b, a[k] as isize, b[k] as isize));
                        sblk[ia * nb + ib] += sq * sv[0] * sv[1] * sv[2];
                        tblk[ia * nb + ib] += sq
                            * (tv[0] * sv[1] * sv[2]
                                + sv[0] * tv[1] * sv[2]
                                + sv[0] * sv[1] * tv[2]);
                    }
                }
                let pref = 2.0 * PI / pp.p * pp.cc;
                for atom in atoms {
                    let c = atom.xyz;
                    r.compute(l, pp.p, pp.pc[0] - c[0], pp.pc[1] - c[1], pp.pc[2] - c[2]);
                    let f = -(atom.z as f64) * pref;
                    for (ia, a) in ca.iter().enumerate() {
                        for (ib, b) in cb.iter().enumerate() {
                            let (ox, oy, oz) = (
                                pp.e[0].offset(a[0], b[0]),
                                pp.e[1].offset(a[1], b[1]),
                                pp.e[2].offset(a[2], b[2]),
                            );
                            let mut acc = 0.0;
                            for tt in 0..=a[0] + b[0] {
                                let ex = pp.e[0].data[ox + tt];
                                for u in 0..=a[1] + b[1] {
                                    let exy = ex * pp.e[1].data[oy + u];
                                    for w in 0..=a[2] + b[2] {
                                        acc += exy * pp.e[2].data[oz + w] * r.at(tt, u, w);
                                    }
                                }
                            }
                            vblk[ia * nb + ib] += f * acc;
                        }
                    }
                }
            }
            scatter(&mut s, &sblk, sa, sb, n);
            scatter(&mut t, &tblk, sa, sb, n);
            scatter(&mut v, &vblk, sa, sb, n);
        }
    }
    OneElectron { s, t, v }
}

/// Dipole integrals <a|r|b> with the origin at 0, as three n x n matrices.
pub fn dipole_matrices(basis: &Basis) -> [Vec<f64>; 3] {
    let n = basis.nbf;
    let mut m = [vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]];
    for si in 0..basis.shells.len() {
        for sj in 0..=si {
            let sa = &basis.shells[si];
            let sb = &basis.shells[sj];
            let ca = cartesian_components(sa.l);
            let cb = cartesian_components(sb.l);
            let nb = cb.len();
            let mut blks = [
                vec![0.0; ca.len() * nb],
                vec![0.0; ca.len() * nb],
                vec![0.0; ca.len() * nb],
            ];
            for pp in prim_pairs(sa, sb, 0, 0) {
                let sq = (PI / pp.p).powf(1.5) * pp.cc;
                for (ia, a) in ca.iter().enumerate() {
                    for (ib, b) in cb.iter().enumerate() {
                        let s: [f64; 3] =
                            [0, 1, 2].map(|k| pp.e[k].get(a[k] as isize, b[k] as isize, 0));
                        for k in 0..3 {
                            let e1 = pp.e[k].get(a[k] as isize, b[k] as isize, 1);
                            let mk = e1 + pp.pc[k] * s[k];
                            blks[k][ia * nb + ib] += sq * mk * s[(k + 1) % 3] * s[(k + 2) % 3];
                        }
                    }
                }
            }
            for k in 0..3 {
                scatter(&mut m[k], &blks[k], sa, sb, n);
            }
        }
    }
    m
}

/// Gradient of the one-electron terms: Σ P dH/dX - Σ W dS/dX (P total density, W energy-weighted
/// density), returned as a 3 * natom vector.
pub fn one_electron_gradient(basis: &Basis, atoms: &[Atom], p: &[f64], w: &[f64]) -> Vec<f64> {
    let n = basis.nbf;
    let mut grad = vec![0.0; atoms.len() * 3];
    let mut r = HermiteR::new();
    for si in 0..basis.shells.len() {
        for sj in 0..=si {
            let sa = &basis.shells[si];
            let sb = &basis.shells[sj];
            let ca = cartesian_components(sa.l);
            let cb = cartesian_components(sb.l);
            let nb = cb.len();
            let fa = component_norms(sa.l);
            let fb = component_norms(sb.l);
            let sym = if si == sj { 1.0 } else { 2.0 };
            let mut pb = vec![0.0; ca.len() * nb];
            let mut wb = vec![0.0; ca.len() * nb];
            for ia in 0..ca.len() {
                for ib in 0..nb {
                    let f = fa[ia] * fb[ib] * sym;
                    let idx = (sa.offset + ia) * n + sb.offset + ib;
                    pb[ia * nb + ib] = p[idx] * f;
                    wb[ia * nb + ib] = w[idx] * f;
                }
            }
            let l = sa.l + sb.l + 1;
            let mut ga = [0.0; 3];
            let mut gb = [0.0; 3];
            let mut gc = vec![0.0; atoms.len() * 3];
            for pp in prim_pairs(sa, sb, 1, 2) {
                let sq = (PI / pp.p).powf(1.5) * pp.cc;
                for (ia, a) in ca.iter().enumerate() {
                    for (ib, b) in cb.iter().enumerate() {
                        let pw = pb[ia * nb + ib];
                        let ww = wb[ia * nb + ib];
                        if pw == 0.0 && ww == 0.0 {
                            continue;
                        }
                        let mut s = [0.0; 3];
                        let mut t = [0.0; 3];
                        let mut ds = [0.0; 3];
                        let mut dt = [0.0; 3];
                        for k in 0..3 {
                            let (i, j) = (a[k] as isize, b[k] as isize);
                            let e = &pp.e[k];
                            s[k] = s1(e, i, j);
                            t[k] = t1(e, pp.b, i, j);
                            ds[k] = 2.0 * pp.a * s1(e, i + 1, j) - i as f64 * s1(e, i - 1, j);
                            dt[k] = 2.0 * pp.a * t1(e, pp.b, i + 1, j)
                                - i as f64 * t1(e, pp.b, i - 1, j);
                        }
                        for k in 0..3 {
                            let (k1, k2) = ((k + 1) % 3, (k + 2) % 3);
                            let d_s = ds[k] * s[k1] * s[k2];
                            let d_t =
                                dt[k] * s[k1] * s[k2] + ds[k] * (t[k1] * s[k2] + s[k1] * t[k2]);
                            let g = sq * (pw * d_t - ww * d_s);
                            ga[k] += g;
                            gb[k] -= g;
                        }
                    }
                }
                let pref = 2.0 * PI / pp.p * pp.cc;
                for (c, atom) in atoms.iter().enumerate() {
                    let cpos = atom.xyz;
                    r.compute(
                        l,
                        pp.p,
                        pp.pc[0] - cpos[0],
                        pp.pc[1] - cpos[1],
                        pp.pc[2] - cpos[2],
                    );
                    let f = -(atom.z as f64) * pref;
                    for (ia, a) in ca.iter().enumerate() {
                        for (ib, b) in cb.iter().enumerate() {
                            let pw = pb[ia * nb + ib];
                            if pw == 0.0 {
                                continue;
                            }
                            let mut da = [0.0; 3];
                            let mut db = [0.0; 3];
                            let (ex_, ey_, ez_) = (&pp.e[0], &pp.e[1], &pp.e[2]);
                            for tt in 0..=a[0] + b[0] + 1 {
                                let ex = ex_.get(a[0] as isize, b[0] as isize, tt);
                                let dax = ex_.d_a(a[0], b[0], tt, pp.a);
                                let dbx = ex_.d_b(a[0], b[0], tt, pp.b);
                                for u in 0..=a[1] + b[1] + 1 {
                                    let ey = ey_.get(a[1] as isize, b[1] as isize, u);
                                    let day = ey_.d_a(a[1], b[1], u, pp.a);
                                    let dby = ey_.d_b(a[1], b[1], u, pp.b);
                                    for vv in 0..=a[2] + b[2] + 1 {
                                        let ez = ez_.get(a[2] as isize, b[2] as isize, vv);
                                        let daz = ez_.d_a(a[2], b[2], vv, pp.a);
                                        let dbz = ez_.d_b(a[2], b[2], vv, pp.b);
                                        let rr = r.at(tt, u, vv);
                                        da[0] += dax * ey * ez * rr;
                                        da[1] += ex * day * ez * rr;
                                        da[2] += ex * ey * daz * rr;
                                        db[0] += dbx * ey * ez * rr;
                                        db[1] += ex * dby * ez * rr;
                                        db[2] += ex * ey * dbz * rr;
                                    }
                                }
                            }
                            let g = f * pw;
                            for k in 0..3 {
                                ga[k] += g * da[k];
                                gb[k] += g * db[k];
                                gc[3 * c + k] -= g * (da[k] + db[k]);
                            }
                        }
                    }
                }
            }
            for k in 0..3 {
                grad[3 * sa.atom + k] += ga[k];
                grad[3 * sb.atom + k] += gb[k];
            }
            for (g, d) in grad.iter_mut().zip(gc.iter()) {
                *g += d;
            }
        }
    }
    grad
}
