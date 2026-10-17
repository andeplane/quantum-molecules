//! Molecules, nuclear repulsion and XYZ parsing.

use crate::elements::{element, element_by_z, BOHR_PER_ANGSTROM};
use crate::Error;

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub z: u32,
    /// Position in bohr.
    pub xyz: [f64; 3],
}

impl Atom {
    pub fn new(symbol: &str, xyz: [f64; 3]) -> Result<Self, Error> {
        let e = element(symbol)
            .ok_or_else(|| Error::new(format!("Unsupported element: {symbol} (H to Ne only)")))?;
        Ok(Atom { z: e.z, xyz })
    }

    pub fn symbol(&self) -> &'static str {
        element_by_z(self.z).map(|e| e.symbol).unwrap_or("?")
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub charge: i32,
    pub multiplicity: u32,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>, charge: i32, multiplicity: u32) -> Self {
        Molecule {
            atoms,
            charge,
            multiplicity,
        }
    }

    pub fn electron_count(&self) -> i64 {
        self.atoms.iter().map(|a| a.z as i64).sum::<i64>() - self.charge as i64
    }

    pub fn nuclear_repulsion(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.atoms.len() {
            for j in 0..i {
                e += (self.atoms[i].z * self.atoms[j].z) as f64
                    / dist(&self.atoms[i].xyz, &self.atoms[j].xyz);
            }
        }
        e
    }

    pub fn nuclear_repulsion_gradient(&self) -> Vec<f64> {
        let n = self.atoms.len();
        let mut g = vec![0.0; 3 * n];
        for i in 0..n {
            for j in 0..i {
                let a = &self.atoms[i];
                let b = &self.atoms[j];
                let r = dist(&a.xyz, &b.xyz);
                let f = -((a.z * b.z) as f64) / (r * r * r);
                for k in 0..3 {
                    let d = a.xyz[k] - b.xyz[k];
                    g[3 * i + k] += f * d;
                    g[3 * j + k] -= f * d;
                }
            }
        }
        g
    }

    /// Flattened coordinates [x0, y0, z0, x1, ...] in bohr.
    pub fn coordinates(&self) -> Vec<f64> {
        self.atoms.iter().flat_map(|a| a.xyz).collect()
    }

    pub fn with_coordinates(&self, x: &[f64]) -> Molecule {
        let mut m = self.clone();
        for (i, a) in m.atoms.iter_mut().enumerate() {
            a.xyz = [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
        }
        m
    }
}

pub fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Parse an XYZ block in Angstrom. Count and comment lines are skipped automatically.
pub fn parse_xyz(text: &str) -> Result<Vec<Atom>, Error> {
    let mut atoms = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let coords: Vec<f64> = match parts[1..4].iter().map(|s| s.parse::<f64>()).collect() {
            Ok(c) => c,
            Err(_) => continue,
        };
        if !parts[0].chars().all(|c| c.is_ascii_alphabetic()) || parts[0].len() > 2 {
            continue;
        }
        atoms.push(Atom::new(
            parts[0],
            [
                coords[0] * BOHR_PER_ANGSTROM,
                coords[1] * BOHR_PER_ANGSTROM,
                coords[2] * BOHR_PER_ANGSTROM,
            ],
        )?);
    }
    if atoms.is_empty() {
        return Err(Error::new("No atoms found in XYZ input"));
    }
    Ok(atoms)
}
