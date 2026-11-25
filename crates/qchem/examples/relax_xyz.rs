//! Relax an XYZ file (Angstrom) at RHF/STO-3G and print the result; used for the larger presets.
//!
//! cargo run --release --example relax_xyz -- molecule.xyz > relaxed.xyz
use qchem::molecule::{parse_xyz, Molecule};
use qchem::optimize::{Bfgs, OptStatus, OptimizerOptions};
use qchem::{run_scf, scf_gradient, BasisName, Method, ScfOptions};
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let atoms = parse_xyz(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let mut mol = Molecule::new(atoms, 0, 1);
    let opts = ScfOptions::new(BasisName::Sto3g, Method::HartreeFock);
    let mut bfgs = Bfgs::new(3 * mol.atoms.len(), OptimizerOptions::default());
    let mut guess = None;
    loop {
        let res = run_scf(&mol, &opts, guess.as_ref()).unwrap();
        let grad = scf_gradient(&mol, &res);
        guess = Some(res.guess());
        match bfgs.step(&mol.coordinates(), res.energy, &grad) {
            OptStatus::Continue(x) => mol = mol.with_coordinates(&x),
            OptStatus::Converged => {
                eprintln!("{path}: {} steps E={:.8}", bfgs.steps, res.energy);
                break;
            }
        }
        if bfgs.steps > 200 {
            eprintln!("{path}: not converged");
            break;
        }
    }
    let b = qchem::elements::BOHR_PER_ANGSTROM;
    for a in &mol.atoms {
        println!(
            "{} {:.4} {:.4} {:.4}",
            a.symbol(),
            a.xyz[0] / b,
            a.xyz[1] / b,
            a.xyz[2] / b
        );
    }
}
