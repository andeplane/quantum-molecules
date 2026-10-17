//! Becke-partitioned molecular integration grids and basis-function evaluation on points.

use crate::basis::{cartesian_components, component_norms, Basis};
use crate::molecule::Atom;
use crate::Error;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GridLevel {
    Coarse,
    Medium,
    Fine,
    Ultrafine,
}

impl GridLevel {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "coarse" => Ok(GridLevel::Coarse),
            "medium" => Ok(GridLevel::Medium),
            "fine" => Ok(GridLevel::Fine),
            "ultrafine" => Ok(GridLevel::Ultrafine),
            _ => Err(Error::new(format!("Unknown grid level: {s}"))),
        }
    }

    /// (radial points for heavy atoms, Lebedev order)
    fn sizes(self) -> (usize, usize) {
        match self {
            GridLevel::Coarse => (35, 110),
            GridLevel::Medium => (50, 194),
            GridLevel::Fine => (75, 302),
            GridLevel::Ultrafine => (120, 590),
        }
    }
}

fn lebedev(order: usize) -> &'static [f64] {
    static DATA: OnceLock<HashMap<String, Vec<f64>>> = OnceLock::new();
    &DATA.get_or_init(|| {
        serde_json::from_str(include_str!("../data/lebedev.json")).expect("valid Lebedev JSON")
    })[&order.to_string()]
}

pub struct MolecularGrid {
    /// x, y, z per point.
    pub coords: Vec<f64>,
    pub weights: Vec<f64>,
    pub npts: usize,
}

/// Mura-Knowles radial quadrature times Lebedev angular grids, Becke fuzzy-cell partitioning.
pub fn molecular_grid(atoms: &[Atom], level: GridLevel) -> MolecularGrid {
    let (radial, angular) = level.sizes();
    let ang = lebedev(angular);
    let nang = ang.len() / 4;
    let natom = atoms.len();
    let mut inv_dist = vec![0.0; natom * natom];
    for i in 0..natom {
        for j in 0..natom {
            if i != j {
                inv_dist[i * natom + j] = 1.0 / crate::molecule::dist(&atoms[i].xyz, &atoms[j].xyz);
            }
        }
    }
    let mut coords = Vec::new();
    let mut weights = Vec::new();
    let mut p = vec![0.0; natom];
    let mut rd = vec![0.0; natom];
    for (ia, atom) in atoms.iter().enumerate() {
        let c = atom.xyz;
        let nr = if atom.z <= 2 {
            (radial as f64 * 0.8).round() as usize
        } else {
            radial
        };
        let alpha_mk = if atom.z == 3 || atom.z == 4 { 7.0 } else { 5.0 };
        for ir in 0..nr {
            let x = (ir as f64 + 0.5) / nr as f64;
            let x3 = x * x * x;
            let r = -alpha_mk * (1.0 - x3).ln();
            let drdx = 3.0 * alpha_mk * x * x / (1.0 - x3);
            let wr = drdx * r * r / nr as f64;
            for k in 0..nang {
                let px = c[0] + r * ang[4 * k];
                let py = c[1] + r * ang[4 * k + 1];
                let pz = c[2] + r * ang[4 * k + 2];
                let w0 = wr * ang[4 * k + 3] * 4.0 * std::f64::consts::PI;
                for i in 0..natom {
                    let a = atoms[i].xyz;
                    rd[i] =
                        ((px - a[0]).powi(2) + (py - a[1]).powi(2) + (pz - a[2]).powi(2)).sqrt();
                    p[i] = 1.0;
                }
                for i in 0..natom {
                    for j in 0..i {
                        let mut mu = (rd[i] - rd[j]) * inv_dist[i * natom + j];
                        for _ in 0..3 {
                            mu = 1.5 * mu - 0.5 * mu * mu * mu;
                        }
                        let s = 0.5 * (1.0 - mu);
                        p[i] *= s;
                        p[j] *= 1.0 - s;
                    }
                }
                let sum: f64 = p.iter().sum();
                let w = w0 * p[ia] / sum;
                if w < 1e-15 {
                    continue;
                }
                coords.extend_from_slice(&[px, py, pz]);
                weights.push(w);
            }
        }
    }
    let npts = weights.len();
    MolecularGrid {
        coords,
        weights,
        npts,
    }
}

/// Basis function values and derivatives on a set of points, layout [g * nbf + mu].
pub struct AoValues {
    pub npts: usize,
    pub nbf: usize,
    pub phi: Vec<f64>,
    /// dx, dy, dz (empty unless deriv >= 1).
    pub d1: [Vec<f64>; 3],
    /// xx, xy, xz, yy, yz, zz (empty unless deriv >= 2).
    pub d2: [Vec<f64>; 6],
}

#[inline(always)]
fn ipow(x: f64, n: usize) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => x * x * x,
        _ => x.powi(n as i32),
    }
}

/// Evaluate basis functions (and up to second derivatives) at `coords` (x, y, z per point).
pub fn eval_ao(basis: &Basis, coords: &[f64], deriv: usize) -> AoValues {
    let npts = coords.len() / 3;
    let n = basis.nbf;
    let alloc = |on: bool| if on { vec![0.0; npts * n] } else { Vec::new() };
    let mut phi = vec![0.0; npts * n];
    let mut d1 = [alloc(deriv >= 1), alloc(deriv >= 1), alloc(deriv >= 1)];
    let mut d2 = [
        alloc(deriv >= 2),
        alloc(deriv >= 2),
        alloc(deriv >= 2),
        alloc(deriv >= 2),
        alloc(deriv >= 2),
        alloc(deriv >= 2),
    ];
    for sh in &basis.shells {
        let comps = cartesian_components(sh.l);
        let norms = component_norms(sh.l);
        let amin = sh.exps.iter().cloned().fold(f64::INFINITY, f64::min);
        for g in 0..npts {
            let x = coords[3 * g] - sh.center[0];
            let y = coords[3 * g + 1] - sh.center[1];
            let z = coords[3 * g + 2] - sh.center[2];
            let r2 = x * x + y * y + z * z;
            if amin * r2 > 46.0 {
                continue;
            }
            let (mut s0, mut s1, mut s2) = (0.0, 0.0, 0.0);
            for (k, &a) in sh.exps.iter().enumerate() {
                let e = sh.coefs[k] * (-a * r2).exp();
                s0 += e;
                s1 += a * e;
                s2 += a * a * e;
            }
            let row = g * n + sh.offset;
            for (c, &[l, m, k]) in comps.iter().enumerate() {
                let f = norms[c];
                let (xl, ym, zk) = (ipow(x, l), ipow(y, m), ipow(z, k));
                phi[row + c] = f * xl * ym * zk * s0;
                if deriv >= 1 {
                    let xl1 = if l > 0 {
                        l as f64 * ipow(x, l - 1)
                    } else {
                        0.0
                    };
                    let ym1 = if m > 0 {
                        m as f64 * ipow(y, m - 1)
                    } else {
                        0.0
                    };
                    let zk1 = if k > 0 {
                        k as f64 * ipow(z, k - 1)
                    } else {
                        0.0
                    };
                    let (xp, yp, zp) = (xl * x, ym * y, zk * z);
                    d1[0][row + c] = f * (xl1 * s0 - 2.0 * xp * s1) * ym * zk;
                    d1[1][row + c] = f * (ym1 * s0 - 2.0 * yp * s1) * xl * zk;
                    d1[2][row + c] = f * (zk1 * s0 - 2.0 * zp * s1) * xl * ym;
                    if deriv >= 2 {
                        let xl2 = if l > 1 {
                            (l * (l - 1)) as f64 * ipow(x, l - 2)
                        } else {
                            0.0
                        };
                        let ym2 = if m > 1 {
                            (m * (m - 1)) as f64 * ipow(y, m - 2)
                        } else {
                            0.0
                        };
                        let zk2 = if k > 1 {
                            (k * (k - 1)) as f64 * ipow(z, k - 2)
                        } else {
                            0.0
                        };
                        let xx = (xl2 * s0 - 2.0 * (2 * l + 1) as f64 * xl * s1
                            + 4.0 * xp * x * s2)
                            * ym
                            * zk;
                        let yy = (ym2 * s0 - 2.0 * (2 * m + 1) as f64 * ym * s1
                            + 4.0 * yp * y * s2)
                            * xl
                            * zk;
                        let zz = (zk2 * s0 - 2.0 * (2 * k + 1) as f64 * zk * s1
                            + 4.0 * zp * z * s2)
                            * xl
                            * ym;
                        let xy = (xl1 * ym1 * s0 - 2.0 * (xl1 * yp + xp * ym1) * s1
                            + 4.0 * xp * yp * s2)
                            * zk;
                        let xz = (xl1 * zk1 * s0 - 2.0 * (xl1 * zp + xp * zk1) * s1
                            + 4.0 * xp * zp * s2)
                            * ym;
                        let yz = (ym1 * zk1 * s0 - 2.0 * (ym1 * zp + yp * zk1) * s1
                            + 4.0 * yp * zp * s2)
                            * xl;
                        for (q, v) in [xx, xy, xz, yy, yz, zz].into_iter().enumerate() {
                            d2[q][row + c] = f * v;
                        }
                    }
                }
            }
        }
    }
    AoValues {
        npts,
        nbf: n,
        phi,
        d1,
        d2,
    }
}
