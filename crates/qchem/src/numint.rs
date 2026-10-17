//! Numerical integration of the exchange-correlation energy, potential and nuclear gradient.

use crate::basis::Basis;
use crate::grid::{eval_ao, MolecularGrid};
use crate::xc::{eval_point, Functional};

const BLOCK: usize = 128;

/// AO values restricted to the basis functions that matter for one block of points.
struct Block {
    start: usize,
    npts: usize,
    idx: Vec<usize>,
    phi: Vec<f64>,
    d1: [Vec<f64>; 3],
    d2: [Vec<f64>; 6],
}

fn blocks<'a>(
    basis: &'a Basis,
    grid: &'a MolecularGrid,
    deriv: usize,
) -> impl Iterator<Item = Block> + 'a {
    let n = basis.nbf;
    (0..grid.npts).step_by(BLOCK).filter_map(move |start| {
        let npts = BLOCK.min(grid.npts - start);
        let ao = eval_ao(basis, &grid.coords[3 * start..3 * (start + npts)], deriv);
        let idx: Vec<usize> = (0..n)
            .filter(|&mu| (0..npts).any(|g| ao.phi[g * n + mu].abs() > 1e-11))
            .collect();
        if idx.is_empty() {
            return None;
        }
        let k = idx.len();
        let pick = |src: &Vec<f64>| -> Vec<f64> {
            if src.is_empty() {
                return Vec::new();
            }
            let mut out = vec![0.0; npts * k];
            for g in 0..npts {
                for (i, &mu) in idx.iter().enumerate() {
                    out[g * k + i] = src[g * n + mu];
                }
            }
            out
        };
        Some(Block {
            start,
            npts,
            phi: pick(&ao.phi),
            d1: [pick(&ao.d1[0]), pick(&ao.d1[1]), pick(&ao.d1[2])],
            d2: [
                pick(&ao.d2[0]),
                pick(&ao.d2[1]),
                pick(&ao.d2[2]),
                pick(&ao.d2[3]),
                pick(&ao.d2[4]),
                pick(&ao.d2[5]),
            ],
            idx,
        })
    })
}

/// T[g, i] = Σ_j A[g, j] D[idx_j, idx_i]
fn contract(a: &[f64], npts: usize, idx: &[usize], d: &[f64], n: usize) -> Vec<f64> {
    let k = idx.len();
    let mut ds = vec![0.0; k * k];
    for i in 0..k {
        for j in 0..k {
            ds[i * k + j] = d[idx[i] * n + idx[j]];
        }
    }
    let mut t = vec![0.0; npts * k];
    for g in 0..npts {
        for j in 0..k {
            let av = a[g * k + j];
            if av == 0.0 {
                continue;
            }
            let row = &ds[j * k..j * k + k];
            let out = &mut t[g * k..g * k + k];
            for i in 0..k {
                out[i] += av * row[i];
            }
        }
    }
    t
}

pub struct XcResult {
    pub exc: f64,
    pub v: Vec<f64>,
    /// Electrons integrated on the grid (a quality check).
    pub nelec: f64,
}

/// XC energy and Kohn-Sham potential matrix for a closed-shell total density `d`.
pub fn xc_matrix(basis: &Basis, grid: &MolecularGrid, d: &[f64], xc: Functional) -> XcResult {
    let n = basis.nbf;
    let gga = xc.is_gga();
    let mut v = vec![0.0; n * n];
    let mut exc = 0.0;
    let mut nelec = 0.0;
    for blk in blocks(basis, grid, if gga { 1 } else { 0 }) {
        let (npts, k) = (blk.npts, blk.idx.len());
        let t = contract(&blk.phi, npts, &blk.idx, d, n);
        let mut x = vec![0.0; npts * k];
        for g in 0..npts {
            let mut rho = 0.0;
            let mut gr = [0.0; 3];
            for i in 0..k {
                let tv = t[g * k + i];
                rho += blk.phi[g * k + i] * tv;
                if gga {
                    for c in 0..3 {
                        gr[c] += blk.d1[c][g * k + i] * tv;
                    }
                }
            }
            for c in gr.iter_mut() {
                *c *= 2.0;
            }
            let sigma = gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2];
            let (e, vrho, vsigma) = eval_point(xc, rho, sigma);
            let w = grid.weights[blk.start + g];
            exc += w * e;
            nelec += w * rho;
            // X[g, i] = w (0.5 vrho phi + 2 vsigma grad rho . grad phi)
            let a = 0.5 * w * vrho;
            let s = 2.0 * w * vsigma;
            for i in 0..k {
                let mut val = a * blk.phi[g * k + i];
                if gga {
                    val += s
                        * (gr[0] * blk.d1[0][g * k + i]
                            + gr[1] * blk.d1[1][g * k + i]
                            + gr[2] * blk.d1[2][g * k + i]);
                }
                x[g * k + i] = val;
            }
        }
        // V += X^T phi + phi^T X
        for i in 0..k {
            for j in 0..=i {
                let mut s = 0.0;
                for g in 0..npts {
                    s += x[g * k + i] * blk.phi[g * k + j] + blk.phi[g * k + i] * x[g * k + j];
                }
                let (mi, mj) = (blk.idx[i], blk.idx[j]);
                v[mi * n + mj] += s;
                if i != j {
                    v[mj * n + mi] += s;
                }
            }
        }
    }
    XcResult { exc, v, nelec }
}

/// XC contribution to the nuclear gradient without grid-weight response (PySCF's default):
/// -2 Σ_{μ∈A,ν} D_μν ∫ [vrho ∂φμ φν + 2 vsigma ∇ρ·(∇∂φμ φν + ∂φμ ∇φν)].
pub fn xc_gradient(
    basis: &Basis,
    grid: &MolecularGrid,
    d: &[f64],
    xc: Functional,
    natom: usize,
) -> Vec<f64> {
    let n = basis.nbf;
    let gga = xc.is_gga();
    let mut grad = vec![0.0; natom * 3];
    const H2: [[usize; 3]; 3] = [[0, 1, 2], [1, 3, 4], [2, 4, 5]];
    for blk in blocks(basis, grid, if gga { 2 } else { 1 }) {
        let (npts, k) = (blk.npts, blk.idx.len());
        let t = contract(&blk.phi, npts, &blk.idx, d, n);
        let mut vrho = vec![0.0; npts];
        let mut vk = [vec![0.0; npts], vec![0.0; npts], vec![0.0; npts]];
        for g in 0..npts {
            let mut rho = 0.0;
            let mut gr = [0.0; 3];
            for i in 0..k {
                let tv = t[g * k + i];
                rho += blk.phi[g * k + i] * tv;
                for c in 0..3 {
                    gr[c] += 2.0 * blk.d1[c][g * k + i] * tv;
                }
            }
            let sigma = gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2];
            let (_, vr, vs) = eval_point(xc, rho, sigma);
            let w = grid.weights[blk.start + g];
            vrho[g] = w * vr;
            for c in 0..3 {
                vk[c][g] = 2.0 * vs * w * gr[c];
            }
        }
        // U = psi D with psi[g, ν] = Σ_k v_k ∂_k φ_ν
        let u = if gga {
            let mut psi = vec![0.0; npts * k];
            for g in 0..npts {
                for i in 0..k {
                    psi[g * k + i] = vk[0][g] * blk.d1[0][g * k + i]
                        + vk[1][g] * blk.d1[1][g * k + i]
                        + vk[2][g] * blk.d1[2][g * k + i];
                }
            }
            contract(&psi, npts, &blk.idx, d, n)
        } else {
            Vec::new()
        };
        for i in 0..k {
            let atom = basis.fn_atom[blk.idx[i]];
            for c in 0..3 {
                let mut s = 0.0;
                for g in 0..npts {
                    let gi = g * k + i;
                    let mut val = vrho[g] * blk.d1[c][gi] * t[gi];
                    if gga {
                        let h = H2[c];
                        val += (vk[0][g] * blk.d2[h[0]][gi]
                            + vk[1][g] * blk.d2[h[1]][gi]
                            + vk[2][g] * blk.d2[h[2]][gi])
                            * t[gi]
                            + blk.d1[c][gi] * u[gi];
                    }
                    s += val;
                }
                grad[3 * atom + c] -= 2.0 * s;
            }
        }
    }
    grad
}
