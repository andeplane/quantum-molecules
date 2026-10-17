//! Contracted cartesian Gaussian basis sets (STO-3G, 6-31G, 6-31G*).

use crate::molecule::Atom;
use crate::Error;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BasisName {
    Sto3g,
    B631g,
    B631gStar,
}

impl BasisName {
    pub fn key(self) -> &'static str {
        match self {
            BasisName::Sto3g => "sto-3g",
            BasisName::B631g => "6-31g",
            BasisName::B631gStar => "6-31g*",
        }
    }

    pub fn parse(s: &str) -> Result<Self, Error> {
        match s.to_ascii_lowercase().as_str() {
            "sto-3g" => Ok(BasisName::Sto3g),
            "6-31g" => Ok(BasisName::B631g),
            "6-31g*" | "6-31g(d)" => Ok(BasisName::B631gStar),
            _ => Err(Error::new(format!("Unknown basis set: {s}"))),
        }
    }
}

#[derive(Deserialize)]
struct RawShell {
    l: usize,
    exps: Vec<f64>,
    coefs: Vec<f64>,
}

type RawData = HashMap<String, HashMap<String, Vec<RawShell>>>;

fn data() -> &'static RawData {
    static DATA: OnceLock<RawData> = OnceLock::new();
    DATA.get_or_init(|| {
        serde_json::from_str(include_str!("../data/basis-data.json")).expect("valid basis JSON")
    })
}

/// A contracted cartesian shell.
#[derive(Clone, Debug)]
pub struct Shell {
    pub atom: usize,
    pub center: [f64; 3],
    pub l: usize,
    pub exps: Vec<f64>,
    /// Contraction coefficients with the radial primitive norm and contraction norm folded in.
    pub coefs: Vec<f64>,
    /// Index of the first basis function of this shell.
    pub offset: usize,
}

#[derive(Clone, Debug)]
pub struct Basis {
    pub name: BasisName,
    pub shells: Vec<Shell>,
    pub nbf: usize,
    /// Atom index for every basis function.
    pub fn_atom: Vec<usize>,
}

pub const fn ncart(l: usize) -> usize {
    (l + 1) * (l + 2) / 2
}

/// Cartesian components in PySCF/libcint order (xx, xy, xz, yy, yz, zz for d).
pub fn cartesian_components(l: usize) -> &'static [[usize; 3]] {
    static TABLE: OnceLock<Vec<Vec<[usize; 3]>>> = OnceLock::new();
    &TABLE.get_or_init(|| {
        (0..=8)
            .map(|l| {
                let mut v = Vec::new();
                for lx in (0..=l).rev() {
                    for ly in (0..=(l - lx)).rev() {
                        v.push([lx, ly, l - lx - ly]);
                    }
                }
                v
            })
            .collect()
    })[l]
}

fn double_factorial(n: i64) -> f64 {
    let mut r = 1.0;
    let mut k = n;
    while k > 1 {
        r *= k as f64;
        k -= 2;
    }
    r
}

/// Factor that turns the (l,0,0)-normalized shell into unit-normalized cartesian components:
/// 1 / sqrt((2lx-1)!! (2ly-1)!! (2lz-1)!!).
pub fn component_norms(l: usize) -> &'static [f64] {
    static TABLE: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    &TABLE.get_or_init(|| {
        (0..=8)
            .map(|l| {
                cartesian_components(l)
                    .iter()
                    .map(|c| {
                        let p: f64 = c
                            .iter()
                            .map(|&k| double_factorial(2 * k as i64 - 1))
                            .product();
                        1.0 / p.sqrt()
                    })
                    .collect()
            })
            .collect()
    })[l]
}

pub fn has_basis(name: BasisName, symbol: &str) -> bool {
    data()
        .get(name.key())
        .is_some_and(|m| m.contains_key(symbol))
}

impl Basis {
    pub fn new(atoms: &[Atom], name: BasisName) -> Result<Basis, Error> {
        let table = &data()[name.key()];
        let mut shells = Vec::new();
        let mut fn_atom = Vec::new();
        let mut offset = 0;
        for (ai, atom) in atoms.iter().enumerate() {
            let sym = atom.symbol();
            let raw = table
                .get(sym)
                .ok_or_else(|| Error::new(format!("No {} basis for {sym}", name.key())))?;
            for rs in raw {
                let l = rs.l;
                let n = rs.exps.len();
                let mut coefs: Vec<f64> = (0..n)
                    .map(|i| {
                        let a = rs.exps[i];
                        rs.coefs[i]
                            * (2.0 * a / std::f64::consts::PI).powf(0.75)
                            * (4.0 * a).powf(l as f64 / 2.0)
                    })
                    .collect();
                let mut s = 0.0;
                for i in 0..n {
                    for j in 0..n {
                        let p = rs.exps[i] + rs.exps[j];
                        s += coefs[i] * coefs[j] * (std::f64::consts::PI / p).powf(1.5)
                            / (2.0 * p).powi(l as i32);
                    }
                }
                let norm = 1.0 / s.sqrt();
                coefs.iter_mut().for_each(|c| *c *= norm);
                shells.push(Shell {
                    atom: ai,
                    center: atom.xyz,
                    l,
                    exps: rs.exps.clone(),
                    coefs,
                    offset,
                });
                fn_atom.extend(std::iter::repeat_n(ai, ncart(l)));
                offset += ncart(l);
            }
        }
        Ok(Basis {
            name,
            shells,
            nbf: offset,
            fn_atom,
        })
    }
}
