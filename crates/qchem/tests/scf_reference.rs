mod common;
use common::*;
use qchem::{run_scf, scf_gradient, GridLevel, ScfOptions};

fn check_case(c: &RefCase) {
    let mol = c.molecule();
    let mut opts = ScfOptions::new(c.basis_name(), c.method_kind());
    opts.grid = GridLevel::Ultrafine;
    opts.conv_energy = 1e-10;
    opts.conv_error = 1e-7;
    let res = run_scf(&mol, &opts, None).unwrap();
    let dft = c.is_dft();
    assert!(res.converged, "{} did not converge", c.label());
    let de = (res.energy - c.energy).abs();
    assert!(
        de < if dft { 2e-5 } else { 1e-8 },
        "{}: energy off by {de}",
        c.label()
    );
    for k in 0..3 {
        let dd = (res.dipole[k] - c.dipole[k]).abs();
        assert!(
            dd < if dft { 1e-4 } else { 1e-6 },
            "{}: dipole off by {dd}",
            c.label()
        );
    }
    if let Some(mo) = &c.mo_energies {
        for (k, e) in mo.iter().enumerate() {
            let d = (res.eps_a[k] - e).abs();
            assert!(
                d < if dft { 1e-4 } else { 1e-6 },
                "{}: orbital {k} off by {d}",
                c.label()
            );
        }
    }
    if let Some(g) = &c.gradient {
        let ours = scf_gradient(&mol, &res);
        let d = max_abs_diff(&ours, &g.concat());
        assert!(
            d < if dft { 1e-4 } else { 1e-6 },
            "{}: gradient off by {d}",
            c.label()
        );
    }
}

macro_rules! reference_test {
    ($name:ident, $mol:expr, $basis:expr, $method:expr, $xc:expr) => {
        #[test]
        fn $name() {
            check_case(&find($mol, $basis, $method, $xc));
        }
    };
}

reference_test!(rhf_water_sto3g, "h2o", "sto-3g", "rhf", None);
reference_test!(rhf_water_631gs, "h2o", "6-31g*", "rhf", None);
reference_test!(
    rhf_distorted_water_631gs,
    "h2o_distorted",
    "6-31g*",
    "rhf",
    None
);
reference_test!(rhf_methane_631gs, "ch4", "6-31g*", "rhf", None);
reference_test!(rhf_ammonia_631g, "nh3", "6-31g", "rhf", None);
reference_test!(rhf_ethylene_sto3g, "c2h4", "sto-3g", "rhf", None);
reference_test!(rhf_nitrogen_631gs, "n2", "6-31g*", "rhf", None);
reference_test!(uhf_hydroxyl_doublet_631g, "oh", "6-31g", "uhf", None);
reference_test!(uhf_methylene_triplet_631gs, "ch2", "6-31g*", "uhf", None);
reference_test!(
    lda_distorted_water_631gs,
    "h2o_distorted",
    "6-31g*",
    "rks",
    Some("lda,vwn")
);
reference_test!(
    pbe_distorted_water_631gs,
    "h2o_distorted",
    "6-31g*",
    "rks",
    Some("pbe")
);
reference_test!(pbe_nitrogen_631g, "n2", "6-31g", "rks", Some("pbe"));
reference_test!(lda_methane_sto3g, "ch4", "sto-3g", "rks", Some("lda,vwn"));

#[test]
fn every_fixture_has_a_test() {
    // Guard against fixtures being added without a matching test above.
    assert_eq!(cases().len(), 13);
}
