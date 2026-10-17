//! Hartree-Fock and PBE on water, then a geometry optimization.
//!
//! cargo run --release --example water

use qchem::molecule::{parse_xyz, Molecule};
use qchem::optimize::{Bfgs, OptStatus, OptimizerOptions};
use qchem::{run_scf, scf_gradient, BasisName, Method, ScfOptions};

fn main() -> Result<(), qchem::Error> {
    let atoms = parse_xyz("O 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692")?;
    let mut mol = Molecule::new(atoms, 0, 1);

    for method in ["hf", "pbe"] {
        let opts = ScfOptions::new(BasisName::B631gStar, Method::parse(method)?);
        let res = run_scf(&mol, &opts, None)?;
        println!(
            "{method}/6-31G*: E = {:.8} Ha in {} iterations",
            res.energy, res.iterations
        );
    }

    let opts = ScfOptions::new(BasisName::Sto3g, Method::HartreeFock);
    let mut bfgs = Bfgs::new(3 * mol.atoms.len(), OptimizerOptions::default());
    let mut guess = None;
    loop {
        let res = run_scf(&mol, &opts, guess.as_ref())?;
        let grad = scf_gradient(&mol, &res);
        guess = Some(res.guess());
        match bfgs.step(&mol.coordinates(), res.energy, &grad) {
            OptStatus::Continue(x) => mol = mol.with_coordinates(&x),
            OptStatus::Converged => {
                println!("relaxed in {} steps: E = {:.8} Ha", bfgs.steps, res.energy);
                break;
            }
        }
    }
    Ok(())
}
