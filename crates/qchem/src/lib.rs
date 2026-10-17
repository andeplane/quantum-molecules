//! `qchem`: a small, dependency-light quantum chemistry engine for molecules made of H to Ne.
//!
//! Hartree-Fock (RHF/UHF) and Kohn-Sham DFT (LDA, PBE) with cartesian Gaussian basis sets
//! (STO-3G, 6-31G, 6-31G*), McMurchie-Davidson integrals, analytic nuclear gradients and a BFGS
//! geometry optimizer. Everything is f64 and single-threaded so it runs unchanged in WebAssembly.

// Index loops mirror the textbook formulas in the integral and linear-algebra kernels.
#![allow(clippy::needless_range_loop, clippy::manual_memcpy)]

pub mod basis;
pub mod boys;
pub mod elements;
pub mod eri;
pub mod gradient;
pub mod grid;
pub mod hermite;
pub mod linalg;
pub mod molecule;
pub mod numint;
pub mod one_electron;
pub mod optimize;
pub mod scf;
pub mod volume;
pub mod xc;

pub use basis::{Basis, BasisName};
pub use gradient::scf_gradient;
pub use grid::GridLevel;
pub use molecule::{Atom, Molecule};
pub use scf::{run_scf, Method, OrbitalGuess, ScfOptions, ScfResult};
pub use xc::Functional;

/// Error type for invalid input or unsupported requests.
#[derive(Clone, Debug, PartialEq)]
pub struct Error(pub String);

impl Error {
    pub fn new(msg: impl Into<String>) -> Self {
        Error(msg.into())
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}
