#![allow(dead_code)]
use qchem::{Atom, BasisName, Method, Molecule};
use serde::Deserialize;

#[derive(Deserialize, Clone)]
pub struct RefAtom {
    pub symbol: String,
    pub xyz: [f64; 3],
}

#[derive(Deserialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct RefCase {
    pub name: String,
    pub basis: String,
    pub method: String,
    pub xc: Option<String>,
    pub charge: i32,
    pub multiplicity: u32,
    pub atoms: Vec<RefAtom>,
    pub nbf: usize,
    pub energy: f64,
    pub dipole: [f64; 3],
    pub mo_energies: Option<Vec<f64>>,
    pub gradient: Option<Vec<[f64; 3]>>,
    #[serde(rename = "S")]
    pub s: Option<Vec<Vec<f64>>>,
    #[serde(rename = "T")]
    pub t: Option<Vec<Vec<f64>>>,
    #[serde(rename = "V")]
    pub v: Option<Vec<Vec<f64>>>,
    pub eri: Option<Vec<f64>>,
    pub eri_index: Option<Vec<usize>>,
}

impl RefCase {
    pub fn molecule(&self) -> Molecule {
        Molecule::new(
            self.atoms
                .iter()
                .map(|a| Atom::new(&a.symbol, a.xyz).unwrap())
                .collect(),
            self.charge,
            self.multiplicity,
        )
    }

    pub fn basis_name(&self) -> BasisName {
        BasisName::parse(&self.basis).unwrap()
    }

    pub fn method_kind(&self) -> Method {
        match (self.method.as_str(), self.xc.as_deref()) {
            ("rks", Some("pbe")) => Method::parse("pbe").unwrap(),
            ("rks", _) => Method::parse("lda").unwrap(),
            _ => Method::HartreeFock,
        }
    }

    pub fn is_dft(&self) -> bool {
        self.method == "rks"
    }

    pub fn label(&self) -> String {
        format!(
            "{} {} {} {}",
            self.name,
            self.basis,
            self.method,
            self.xc.clone().unwrap_or_default()
        )
    }
}

pub fn cases() -> Vec<RefCase> {
    serde_json::from_str(include_str!("../fixtures/pyscf-reference.json")).unwrap()
}

pub fn find(name: &str, basis: &str, method: &str, xc: Option<&str>) -> RefCase {
    cases()
        .into_iter()
        .find(|c| c.name == name && c.basis == basis && c.method == method && c.xc.as_deref() == xc)
        .unwrap_or_else(|| panic!("missing fixture {name} {basis} {method}"))
}

pub fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Water near its RHF/STO-3G minimum, in bohr.
pub fn water() -> Molecule {
    let ang = qchem::elements::BOHR_PER_ANGSTROM;
    Molecule::new(
        vec![
            Atom::new("O", [0.0, 0.0, 0.1173 * ang]).unwrap(),
            Atom::new("H", [0.0, 0.7572 * ang, -0.4692 * ang]).unwrap(),
            Atom::new("H", [0.0, -0.7572 * ang, -0.4692 * ang]).unwrap(),
        ],
        0,
        1,
    )
}
