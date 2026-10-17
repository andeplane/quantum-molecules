//! Analytic nuclear gradients of SCF energies.

use crate::eri::eri_gradient;
use crate::molecule::Molecule;
use crate::numint::xc_gradient;
use crate::one_electron::one_electron_gradient;
use crate::scf::{Method, ScfResult};

/// dE/dR in Hartree/bohr, flattened as [x0, y0, z0, x1, ...].
pub fn scf_gradient(mol: &Molecule, res: &ScfResult) -> Vec<f64> {
    let n = res.nbf;
    let nmo = res.nmo;
    let natom = mol.atoms.len();
    let dt = res.total_density();
    // Energy-weighted density W = Σ_σ Σ_occ ε_i C_μi C_νi
    let mut w = vec![0.0; n * n];
    for (c, eps, nocc) in [
        (&res.ca, &res.eps_a, res.nalpha),
        (&res.cb, &res.eps_b, res.nbeta),
    ] {
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..nocc {
                    s += eps[k] * c[i * nmo + k] * c[j * nmo + k];
                }
                w[i * n + j] += s;
            }
        }
    }
    let kx = res.method.exact_exchange();
    let mut grad = mol.nuclear_repulsion_gradient();
    let g1 = one_electron_gradient(&res.basis, &mol.atoms, &dt, &w);
    let g2 = eri_gradient(&res.basis, natom, &dt, &res.da, &res.db, kx);
    for i in 0..grad.len() {
        grad[i] += g1[i] + g2[i];
    }
    if let (Method::Dft(f), Some(grid)) = (res.method, &res.grid) {
        let g3 = xc_gradient(&res.basis, grid, &dt, f, natom);
        for i in 0..grad.len() {
            grad[i] += g3[i];
        }
    }
    grad
}
