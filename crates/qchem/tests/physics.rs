//! Physical consistency checks that do not need external reference data.
mod common;
use common::*;
use qchem::elements::BOHR_PER_ANGSTROM;
use qchem::grid::molecular_grid;
use qchem::molecule::{dist, parse_xyz};
use qchem::optimize::{Bfgs, OptStatus, OptimizerOptions};
use qchem::volume::{density_at, VolumeBox};
use qchem::{run_scf, scf_gradient, Atom, BasisName, GridLevel, Method, Molecule, ScfOptions};

fn opts(basis: BasisName, method: &str) -> ScfOptions {
    let mut o = ScfOptions::new(basis, Method::parse(method).unwrap());
    o.conv_energy = 1e-11;
    o.conv_error = 1e-8;
    o.grid = GridLevel::Fine;
    o
}

fn energy(mol: &Molecule, o: &ScfOptions) -> f64 {
    let r = run_scf(mol, o, None).unwrap();
    assert!(r.converged);
    r.energy
}

fn distorted_water() -> Molecule {
    let mut m = water();
    m.atoms[0].xyz[0] += 0.07;
    m.atoms[1].xyz[1] += 0.12;
    m.atoms[2].xyz[2] -= 0.09;
    m
}

/// Central finite differences of the energy versus the analytic gradient.
fn check_fd(mol: &Molecule, o: &ScfOptions, tol: f64) {
    let res = run_scf(mol, o, None).unwrap();
    let g = scf_gradient(mol, &res);
    let x = mol.coordinates();
    let h = 1e-3;
    for i in 0..x.len() {
        let mut xp = x.clone();
        xp[i] += h;
        let mut xm = x.clone();
        xm[i] -= h;
        let fd = (energy(&mol.with_coordinates(&xp), o) - energy(&mol.with_coordinates(&xm), o))
            / (2.0 * h);
        assert!(
            (fd - g[i]).abs() < tol,
            "coordinate {i}: analytic {} vs finite difference {fd}",
            g[i]
        );
    }
}

#[test]
fn rhf_gradient_matches_finite_differences() {
    check_fd(&distorted_water(), &opts(BasisName::B631gStar, "hf"), 2e-6);
}

#[test]
fn uhf_gradient_matches_finite_differences() {
    let ang = BOHR_PER_ANGSTROM;
    let oh = Molecule::new(
        vec![
            Atom::new("O", [0.0, 0.0, 0.0]).unwrap(),
            Atom::new("H", [0.1, 0.05, 1.02 * ang]).unwrap(),
        ],
        0,
        2,
    );
    check_fd(&oh, &opts(BasisName::B631g, "hf"), 2e-6);
}

#[test]
fn dft_gradient_matches_finite_differences() {
    // The grid moves with the atoms but the analytic gradient ignores grid-weight derivatives,
    // so agreement is limited by grid quality.
    check_fd(&distorted_water(), &opts(BasisName::Sto3g, "pbe"), 3e-4);
    check_fd(&distorted_water(), &opts(BasisName::Sto3g, "lda"), 3e-4);
}

fn rotate(mol: &Molecule, angle: f64) -> Molecule {
    let (c, s) = (angle.cos(), angle.sin());
    let mut m = mol.clone();
    for a in m.atoms.iter_mut() {
        let [x, y, z] = a.xyz;
        // rotation about (1,1,1)/sqrt(3) is overkill; compose z then x rotations.
        let (x1, y1) = (c * x - s * y, s * x + c * y);
        let (y2, z2) = (c * y1 - s * z, s * y1 + c * z);
        a.xyz = [x1 + 0.37, y2 - 1.1, z2 + 2.3];
    }
    m
}

#[test]
fn energy_is_invariant_to_rotation_and_translation() {
    let mol = distorted_water();
    for method in ["hf", "pbe"] {
        let o = opts(BasisName::B631gStar, method);
        let e0 = energy(&mol, &o);
        let e1 = energy(&rotate(&mol, 0.7), &o);
        // Cartesian d functions are rotation-invariant as a set; DFT grids are not exactly.
        let tol = if method == "hf" { 1e-9 } else { 2e-5 };
        assert!((e0 - e1).abs() < tol, "{method}: {e0} vs {e1}");
    }
}

#[test]
fn gradient_has_no_net_force() {
    let mol = distorted_water();
    for method in ["hf", "lda"] {
        let o = opts(BasisName::B631gStar, method);
        let res = run_scf(&mol, &o, None).unwrap();
        let g = scf_gradient(&mol, &res);
        for k in 0..3 {
            let net: f64 = (0..mol.atoms.len()).map(|a| g[3 * a + k]).sum();
            let tol = if method == "hf" { 1e-9 } else { 1e-4 };
            assert!(net.abs() < tol, "{method}: net force {net}");
        }
    }
}

#[test]
fn warm_start_converges_faster_to_the_same_state() {
    let mol = distorted_water();
    let o = opts(BasisName::B631gStar, "hf");
    let first = run_scf(&mol, &o, None).unwrap();
    let mut moved = mol.clone();
    moved.atoms[1].xyz[1] += 0.02;
    let cold = run_scf(&moved, &o, None).unwrap();
    let warm = run_scf(&moved, &o, Some(&first.guess())).unwrap();
    assert!((cold.energy - warm.energy).abs() < 1e-9);
    assert!(
        warm.iterations < cold.iterations,
        "warm {} cold {}",
        warm.iterations,
        cold.iterations
    );
}

#[test]
fn history_records_every_iteration() {
    let mut o = opts(BasisName::Sto3g, "hf");
    o.record_history = true;
    let r = run_scf(&water(), &o, None).unwrap();
    assert_eq!(r.history.len(), r.iterations);
    let last = r.history.last().unwrap();
    assert!((last.energy - r.energy).abs() < 1e-12);
    assert_eq!(last.ca.len(), r.nbf * r.nmo);
}

#[test]
fn charged_molecules_are_supported() {
    // Hydroxide anion: 10 electrons, closed shell.
    let ang = BOHR_PER_ANGSTROM;
    let oh = Molecule::new(
        vec![
            Atom::new("O", [0.0, 0.0, 0.0]).unwrap(),
            Atom::new("H", [0.0, 0.0, 0.97 * ang]).unwrap(),
        ],
        -1,
        1,
    );
    let r = run_scf(&oh, &opts(BasisName::B631gStar, "hf"), None).unwrap();
    assert!(r.converged);
    assert_eq!(r.nalpha + r.nbeta, 10);
    // An anion's HOMO is bound but shallow; a sanity window.
    assert!(r.eps_a[4] < 0.0 && r.eps_a[4] > -0.3);
}

#[test]
fn invalid_requests_are_rejected() {
    let w = water();
    let o = opts(BasisName::Sto3g, "hf");
    let mut bad = w.clone();
    bad.multiplicity = 2; // 10 electrons cannot be a doublet
    assert!(run_scf(&bad, &o, None).is_err());
    let mut radical = w.clone();
    radical.charge = 1;
    radical.multiplicity = 2;
    assert!(run_scf(&radical, &opts(BasisName::Sto3g, "pbe"), None).is_err());
    assert!(run_scf(&radical, &o, None).unwrap().unrestricted);
    assert!(Atom::new("Cl", [0.0; 3]).is_err());
    assert!(BasisName::parse("cc-pvdz").is_err());
    assert!(Method::parse("b3lyp").is_err());
    assert!(run_scf(&Molecule::new(vec![], 0, 1), &o, None).is_err());
}

#[test]
fn relaxes_water_to_the_known_sto3g_minimum() {
    // RHF/STO-3G water: r(OH) = 0.989 A, angle(HOH) = 100.0 degrees.
    let mut mol = distorted_water();
    let o = opts(BasisName::Sto3g, "hf");
    let mut bfgs = Bfgs::new(
        9,
        OptimizerOptions {
            max_force: 1e-5,
            rms_force: 1e-5,
            ..Default::default()
        },
    );
    let mut guess = None;
    let start = energy(&mol, &o);
    let mut converged = false;
    for _ in 0..80 {
        let r = run_scf(&mol, &o, guess.as_ref()).unwrap();
        let g = scf_gradient(&mol, &r);
        guess = Some(r.guess());
        match bfgs.step(&mol.coordinates(), r.energy, &g) {
            OptStatus::Continue(x) => mol = mol.with_coordinates(&x),
            OptStatus::Converged => {
                converged = true;
                break;
            }
        }
    }
    assert!(converged);
    assert!(energy(&mol, &o) < start);
    let r1 = dist(&mol.atoms[0].xyz, &mol.atoms[1].xyz) / BOHR_PER_ANGSTROM;
    let r2 = dist(&mol.atoms[0].xyz, &mol.atoms[2].xyz) / BOHR_PER_ANGSTROM;
    let v1: Vec<f64> = (0..3)
        .map(|k| mol.atoms[1].xyz[k] - mol.atoms[0].xyz[k])
        .collect();
    let v2: Vec<f64> = (0..3)
        .map(|k| mol.atoms[2].xyz[k] - mol.atoms[0].xyz[k])
        .collect();
    let cos =
        v1.iter().zip(&v2).map(|(a, b)| a * b).sum::<f64>() / (r1 * r2 * BOHR_PER_ANGSTROM.powi(2));
    let angle = cos.acos().to_degrees();
    assert!(
        (r1 - 0.989).abs() < 0.002 && (r2 - 0.989).abs() < 0.002,
        "r = {r1}, {r2}"
    );
    assert!((angle - 100.0).abs() < 0.3, "angle = {angle}");
}

#[test]
fn density_integrates_to_the_electron_count() {
    let mol = distorted_water();
    let r = run_scf(&mol, &opts(BasisName::B631gStar, "hf"), None).unwrap();
    let grid = molecular_grid(&mol.atoms, GridLevel::Fine);
    let rho = density_at(&r, &grid.coords);
    let n: f64 = rho.iter().zip(&grid.weights).map(|(a, b)| a * b).sum();
    assert!((n - 10.0).abs() < 1e-4, "{n}");
}

#[test]
fn volume_box_spans_the_molecule() {
    let mol = water();
    let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.xyz).collect();
    let b = VolumeBox::around(&coords, 4.0, [3, 4, 5]);
    let pts = b.points();
    assert_eq!(pts.len(), 3 * 60);
    assert_eq!(&pts[0..3], &b.min);
    assert_eq!(&pts[pts.len() - 3..], &b.max);
}

#[test]
fn parses_xyz_input() {
    let atoms =
        parse_xyz("3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nh 0 -0.7572 -0.4692\n").unwrap();
    assert_eq!(atoms.len(), 3);
    assert_eq!(atoms[2].symbol(), "H");
    assert!((atoms[0].xyz[2] - 0.1173 * BOHR_PER_ANGSTROM).abs() < 1e-12);
    assert!(parse_xyz("nothing here").is_err());
    assert!(parse_xyz("Xx 0 0 0").is_err());
}

/// <previous_a | new_b> for the orbitals of two results, using the new geometry's overlap.
fn orbital_overlap(prev: &qchem::ScfResult, new: &qchem::ScfResult) -> Vec<f64> {
    use qchem::linalg::{matmul, transpose};
    let s = qchem::one_electron::one_electron_matrices(&new.basis, &[]).s;
    let (n, m) = (new.nbf, new.nmo);
    matmul(
        &transpose(&prev.ca, n, m),
        &matmul(&s, &new.ca, n, n, m),
        m,
        n,
        m,
    )
}

#[test]
fn orbital_phases_follow_a_dragged_atom() {
    // Dragging a hydrogen in small steps, as the app does: every orbital must keep its sign,
    // otherwise the orbital isosurface flips colour from frame to frame.
    let o = opts(BasisName::B631g, "hf");
    let mut mol = water();
    let mut prev = run_scf(&mol, &o, None).unwrap();
    for step in 0..12 {
        mol.atoms[1].xyz[1] += 0.03;
        mol.atoms[1].xyz[2] += if step % 2 == 0 { 0.02 } else { -0.01 };
        let next = run_scf(&mol, &o, Some(&prev.guess())).unwrap();
        let ov = orbital_overlap(&prev, &next);
        for k in 0..next.nalpha + 2 {
            assert!(
                ov[k * next.nmo + k] > 0.9,
                "step {step}, orbital {k}: overlap {}",
                ov[k * next.nmo + k]
            );
        }
        prev = next;
    }
}

#[test]
fn cold_starts_give_identical_orbitals() {
    let o = opts(BasisName::B631gStar, "hf");
    let a = run_scf(&distorted_water(), &o, None).unwrap();
    let b = run_scf(&distorted_water(), &o, None).unwrap();
    assert!(max_abs_diff(&a.ca, &b.ca) < 1e-12);
    // Convention without a reference: the largest coefficient of each orbital is positive.
    for k in 0..a.nmo {
        let col: Vec<f64> = (0..a.nbf).map(|i| a.ca[i * a.nmo + k]).collect();
        let big = col
            .iter()
            .cloned()
            .fold(0.0_f64, |m, v| if v.abs() > m.abs() { v } else { m });
        assert!(big > 0.0, "orbital {k}");
    }
}

#[test]
fn degenerate_orbitals_follow_the_previous_ones() {
    // N2's pi orbitals come in degenerate pairs that any rotation of solves equally well; after
    // stretching the bond they must stay aligned with the previous pair, not spin around.
    let o = opts(BasisName::B631gStar, "hf");
    let ang = BOHR_PER_ANGSTROM;
    let n2 = |r: f64| {
        Molecule::new(
            vec![
                Atom::new("N", [0.0, 0.0, -r / 2.0 * ang]).unwrap(),
                Atom::new("N", [0.0, 0.0, r / 2.0 * ang]).unwrap(),
            ],
            0,
            1,
        )
    };
    let prev = run_scf(&n2(1.098), &o, None).unwrap();
    let next = run_scf(&n2(1.12), &o, Some(&prev.guess())).unwrap();
    let ov = orbital_overlap(&prev, &next);
    let m = next.nmo;
    for k in 0..next.nalpha + 2 {
        assert!(
            ov[k * m + k] > 0.99,
            "orbital {k}: overlap {}",
            ov[k * m + k]
        );
    }
    // The energy and density are untouched by the alignment.
    let cold = run_scf(&n2(1.12), &o, None).unwrap();
    assert!((cold.energy - next.energy).abs() < 1e-9);
    assert!(max_abs_diff(&cold.da, &next.da) < 1e-6);
}
