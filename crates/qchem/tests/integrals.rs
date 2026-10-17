mod common;
use common::*;
use qchem::basis::Basis;
use qchem::eri::{build_jk, compute_eris};
use qchem::one_electron::one_electron_matrices;

#[test]
fn basis_function_counts_match_pyscf() {
    for c in cases() {
        let b = Basis::new(&c.molecule().atoms, c.basis_name()).unwrap();
        assert_eq!(b.nbf, c.nbf, "{}", c.label());
    }
}

#[test]
fn one_electron_matrices_match_pyscf() {
    for basis in ["sto-3g", "6-31g*"] {
        let c = find("h2o", basis, "rhf", None);
        let mol = c.molecule();
        let b = Basis::new(&mol.atoms, c.basis_name()).unwrap();
        let m = one_electron_matrices(&b, &mol.atoms);
        let flat = |x: &Option<Vec<Vec<f64>>>| x.as_ref().unwrap().concat();
        assert!(max_abs_diff(&m.s, &flat(&c.s)) < 1e-10, "S {basis}");
        assert!(max_abs_diff(&m.t, &flat(&c.t)) < 1e-9, "T {basis}");
        assert!(max_abs_diff(&m.v, &flat(&c.v)) < 1e-9, "V {basis}");
    }
}

#[test]
fn basis_functions_are_normalized() {
    let c = find("n2", "6-31g*", "rhf", None);
    let mol = c.molecule();
    let b = Basis::new(&mol.atoms, c.basis_name()).unwrap();
    let s = one_electron_matrices(&b, &mol.atoms).s;
    for i in 0..b.nbf {
        assert!((s[i * b.nbf + i] - 1.0).abs() < 1e-12);
    }
}

#[test]
fn electron_repulsion_integrals_match_pyscf() {
    for basis in ["sto-3g", "6-31g*"] {
        let c = find("h2o", basis, "rhf", None);
        let mol = c.molecule();
        let b = Basis::new(&mol.atoms, c.basis_name()).unwrap();
        let eri = compute_eris(&b);
        let n = b.nbf;
        let mut worst: f64 = 0.0;
        for (k, &flat) in c.eri_index.as_ref().unwrap().iter().enumerate() {
            let (i, j, kk, l) = (
                flat / (n * n * n),
                (flat / (n * n)) % n,
                (flat / n) % n,
                flat % n,
            );
            worst = worst.max((eri.get(i, j, kk, l) - c.eri.as_ref().unwrap()[k]).abs());
        }
        assert!(worst < 1e-10, "{basis}: {worst}");
    }
}

#[test]
fn eri_tensor_has_eightfold_symmetry() {
    let c = find("h2o", "6-31g*", "rhf", None);
    let mol = c.molecule();
    let b = Basis::new(&mol.atoms, c.basis_name()).unwrap();
    let eri = compute_eris(&b);
    // Packed storage returns the same element for all 8 permutations by construction; check the
    // index map itself and positivity of the diagonal (ij|ij) >= 0.
    for (i, j, k, l) in [(3, 1, 7, 2), (10, 4, 4, 9), (18, 0, 5, 5)] {
        let v = eri.get(i, j, k, l);
        for (p, q, r, s) in [(j, i, k, l), (i, j, l, k), (k, l, i, j), (l, k, j, i)] {
            assert_eq!(eri.get(p, q, r, s), v);
        }
    }
    for i in 0..b.nbf {
        for j in 0..b.nbf {
            assert!(eri.get(i, j, i, j) >= -1e-14);
        }
    }
}

#[test]
fn coulomb_and_exchange_match_dense_contraction() {
    let c = find("h2o", "sto-3g", "rhf", None);
    let mol = c.molecule();
    let b = Basis::new(&mol.atoms, c.basis_name()).unwrap();
    let n = b.nbf;
    let eri = compute_eris(&b);
    // A non-symmetric-looking but symmetric test density.
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = ((i * 7 + j * 3) % 5) as f64 * 0.1 - 0.2;
            d[i * n + j] = v;
            d[j * n + i] = v;
        }
    }
    let (j, k) = build_jk(&eri, &d, &[&d]);
    for p in 0..n {
        for q in 0..n {
            let mut jj = 0.0;
            let mut kk = 0.0;
            for r in 0..n {
                for s in 0..n {
                    jj += eri.get(p, q, r, s) * d[r * n + s];
                    kk += eri.get(p, r, q, s) * d[r * n + s];
                }
            }
            assert!((j[p * n + q] - jj).abs() < 1e-12);
            assert!((k[0][p * n + q] - kk).abs() < 1e-12);
        }
    }
}
