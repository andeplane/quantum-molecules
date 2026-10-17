//! Electron density and orbital values on points or regular 3D boxes (CPU reference for the
//! GPU volume renderer, and for exporting cube-like data).

use crate::basis::Basis;
use crate::grid::eval_ao;
use crate::scf::ScfResult;

/// An axis-aligned box sampled at `dims` points per axis (bohr).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub dims: [usize; 3],
}

impl VolumeBox {
    /// Bounding box of the atoms padded by `margin` bohr on every side.
    pub fn around(coords: &[[f64; 3]], margin: f64, dims: [usize; 3]) -> Self {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for c in coords {
            for k in 0..3 {
                min[k] = min[k].min(c[k] - margin);
                max[k] = max[k].max(c[k] + margin);
            }
        }
        VolumeBox { min, max, dims }
    }

    pub fn point(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let idx = [i, j, k];
        let mut p = [0.0; 3];
        for a in 0..3 {
            let t = if self.dims[a] > 1 {
                idx[a] as f64 / (self.dims[a] - 1) as f64
            } else {
                0.5
            };
            p[a] = self.min[a] + t * (self.max[a] - self.min[a]);
        }
        p
    }

    /// All sample points, x fastest.
    pub fn points(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(3 * self.dims.iter().product::<usize>());
        for k in 0..self.dims[2] {
            for j in 0..self.dims[1] {
                for i in 0..self.dims[0] {
                    out.extend_from_slice(&self.point(i, j, k));
                }
            }
        }
        out
    }
}

/// Values of molecular orbital columns `orbs` of `c` (nbf x nmo) at `coords`, layout [orb][point].
pub fn orbitals_at(
    basis: &Basis,
    c: &[f64],
    nmo: usize,
    orbs: &[usize],
    coords: &[f64],
) -> Vec<Vec<f64>> {
    let n = basis.nbf;
    let npts = coords.len() / 3;
    let mut out = vec![vec![0.0; npts]; orbs.len()];
    for chunk in (0..npts).step_by(256) {
        let m = 256.min(npts - chunk);
        let ao = eval_ao(basis, &coords[3 * chunk..3 * (chunk + m)], 0);
        for (o, &k) in orbs.iter().enumerate() {
            for g in 0..m {
                let mut s = 0.0;
                for mu in 0..n {
                    s += ao.phi[g * n + mu] * c[mu * nmo + k];
                }
                out[o][chunk + g] = s;
            }
        }
    }
    out
}

/// Total electron density at `coords`.
pub fn density_at(res: &ScfResult, coords: &[f64]) -> Vec<f64> {
    let npts = coords.len() / 3;
    let mut rho = vec![0.0; npts];
    let spins: [(&[f64], usize); 2] = [(&res.ca, res.nalpha), (&res.cb, res.nbeta)];
    for (c, nocc) in spins {
        let orbs: Vec<usize> = (0..nocc).collect();
        for psi in orbitals_at(&res.basis, c, res.nmo, &orbs, coords) {
            for (r, v) in rho.iter_mut().zip(psi) {
                *r += v * v;
            }
        }
    }
    rho
}
