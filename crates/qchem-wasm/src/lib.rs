//! WebAssembly bindings for `qchem`.
//!
//! JavaScript passes plain objects (deserialized with serde) and receives plain objects back.
//! `Session` keeps the last converged orbitals so nearby geometries start from them.

use qchem::basis::cartesian_components;
use qchem::molecule::Molecule;
use qchem::optimize::{force_stats, Bfgs, OptStatus, OptimizerOptions};
use qchem::{
    run_scf, scf_gradient, Atom, BasisName, GridLevel, Method, OrbitalGuess, ScfOptions, ScfResult,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AtomInput {
    pub symbol: String,
    /// Bohr.
    pub xyz: [f64; 3],
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ComputeInput {
    pub atoms: Vec<AtomInput>,
    #[serde(default)]
    pub charge: i32,
    #[serde(default = "one")]
    pub multiplicity: u32,
    pub basis: String,
    pub method: String,
    #[serde(default)]
    pub grid: Option<String>,
    #[serde(default)]
    pub gradient: bool,
    #[serde(default)]
    pub history: bool,
    /// Start from the previous result's orbitals when the basis matches.
    #[serde(default = "yes")]
    pub warm_start: bool,
}

fn one() -> u32 {
    1
}
fn yes() -> bool {
    true
}

/// Basis set laid out for a GPU shader: per shell [cx, cy, cz, l, nprim, primStart, fnOffset, 0]
/// and per primitive [exponent, coefficient]. Coefficients include all normalization except the
/// per-component factor, which is `componentNorms[l][k]`.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GpuBasis {
    pub shells: Vec<f64>,
    pub prims: Vec<f64>,
    pub component_norms: Vec<Vec<f64>>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct IterationOutput {
    pub energy: f64,
    pub error: f64,
    pub ca: Vec<f64>,
    pub cb: Vec<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ComputeOutput {
    pub energy: f64,
    pub converged: bool,
    pub iterations: usize,
    pub unrestricted: bool,
    pub nbf: usize,
    pub nmo: usize,
    pub nalpha: usize,
    pub nbeta: usize,
    pub eps_a: Vec<f64>,
    pub eps_b: Vec<f64>,
    /// nbf x nmo, row-major.
    pub ca: Vec<f64>,
    pub cb: Vec<f64>,
    /// Atomic units (e bohr).
    pub dipole: [f64; 3],
    pub nuclear_repulsion: f64,
    pub exc: f64,
    pub gradient: Option<Vec<f64>>,
    pub max_force: Option<f64>,
    pub history: Vec<IterationOutput>,
    pub basis: GpuBasis,
    /// Atom positions the result belongs to (bohr), flattened.
    pub coordinates: Vec<f64>,
    pub millis: f64,
}

fn molecule_from(input: &ComputeInput) -> Result<Molecule, JsError> {
    let atoms = input
        .atoms
        .iter()
        .map(|a| Atom::new(&a.symbol, a.xyz))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| JsError::new(&e.0))?;
    Ok(Molecule::new(atoms, input.charge, input.multiplicity))
}

fn gpu_basis(res: &ScfResult) -> GpuBasis {
    let mut shells = Vec::new();
    let mut prims = Vec::new();
    for sh in &res.basis.shells {
        let start = prims.len() / 2;
        for (a, c) in sh.exps.iter().zip(&sh.coefs) {
            prims.push(*a);
            prims.push(*c);
        }
        shells.extend_from_slice(&[
            sh.center[0],
            sh.center[1],
            sh.center[2],
            sh.l as f64,
            sh.exps.len() as f64,
            start as f64,
            sh.offset as f64,
            0.0,
        ]);
    }
    let component_norms = (0..=2)
        .map(|l| qchem::basis::component_norms(l).to_vec())
        .collect::<Vec<_>>();
    debug_assert!(cartesian_components(2).len() == 6);
    GpuBasis {
        shells,
        prims,
        component_norms,
    }
}

fn now() -> f64 {
    js_sys::Date::now()
}

/// Holds the previous result so the next geometry can start from its orbitals.
#[wasm_bindgen]
pub struct Session {
    guess: Option<OrbitalGuess>,
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

#[wasm_bindgen]
impl Session {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Session {
        Session { guess: None }
    }

    /// Forget the stored orbitals (for example after switching molecules).
    pub fn reset(&mut self) {
        self.guess = None;
    }

    /// Run an SCF calculation. See `ComputeInput` for the fields of `input`.
    pub fn compute(&mut self, input: JsValue) -> Result<JsValue, JsError> {
        let t0 = now();
        let input: ComputeInput = serde_wasm_bindgen::from_value(input)?;
        let mol = molecule_from(&input)?;
        let basis = BasisName::parse(&input.basis).map_err(|e| JsError::new(&e.0))?;
        let method = Method::parse(&input.method).map_err(|e| JsError::new(&e.0))?;
        let mut opts = ScfOptions::new(basis, method);
        opts.record_history = input.history;
        if let Some(g) = &input.grid {
            opts.grid = GridLevel::parse(g).map_err(|e| JsError::new(&e.0))?;
        }
        let guess = if input.warm_start {
            self.guess.as_ref()
        } else {
            None
        };
        let mut res = run_scf(&mol, &opts, guess).map_err(|e| JsError::new(&e.0))?;
        if !res.converged && guess.is_some() {
            // A stale guess can trap the SCF; retry from atomic densities.
            res = run_scf(&mol, &opts, None).map_err(|e| JsError::new(&e.0))?;
        }
        if res.converged {
            self.guess = Some(res.guess());
        }
        let gradient = if input.gradient {
            Some(scf_gradient(&mol, &res))
        } else {
            None
        };
        let out = ComputeOutput {
            energy: res.energy,
            converged: res.converged,
            iterations: res.iterations,
            unrestricted: res.unrestricted,
            nbf: res.nbf,
            nmo: res.nmo,
            nalpha: res.nalpha,
            nbeta: res.nbeta,
            eps_a: res.eps_a.clone(),
            eps_b: res.eps_b.clone(),
            ca: res.ca.clone(),
            cb: res.cb.clone(),
            dipole: res.dipole,
            nuclear_repulsion: res.nuclear_repulsion,
            exc: res.exc,
            max_force: gradient.as_ref().map(|g| force_stats(g).0),
            gradient,
            history: res
                .history
                .iter()
                .map(|h| IterationOutput {
                    energy: h.energy,
                    error: h.error,
                    ca: h.ca.clone(),
                    cb: h.cb.clone(),
                })
                .collect(),
            basis: gpu_basis(&res),
            coordinates: mol.coordinates(),
            millis: now() - t0,
        };
        Ok(serde_wasm_bindgen::to_value(&out)?)
    }
}

/// Step-wise BFGS geometry optimizer.
#[wasm_bindgen]
pub struct Optimizer {
    inner: Bfgs,
}

#[wasm_bindgen]
impl Optimizer {
    #[wasm_bindgen(constructor)]
    pub fn new(natom: usize) -> Optimizer {
        Optimizer {
            inner: Bfgs::new(3 * natom, OptimizerOptions::default()),
        }
    }

    /// Report energy and gradient at `coords` (bohr). Returns the next coordinates to evaluate,
    /// or `undefined` once converged.
    pub fn step(&mut self, coords: Vec<f64>, energy: f64, gradient: Vec<f64>) -> Option<Vec<f64>> {
        match self.inner.step(&coords, energy, &gradient) {
            OptStatus::Continue(x) => Some(x),
            OptStatus::Converged => None,
        }
    }

    #[wasm_bindgen(getter)]
    pub fn steps(&self) -> usize {
        self.inner.steps
    }
}

/// Which elements a basis set covers, for the UI.
#[wasm_bindgen(js_name = hasBasis)]
pub fn has_basis(basis: &str, symbol: &str) -> bool {
    BasisName::parse(basis)
        .map(|b| qchem::basis::has_basis(b, symbol))
        .unwrap_or(false)
}
