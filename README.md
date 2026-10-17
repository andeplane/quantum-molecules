# Quantum Molecules

Real quantum chemistry for small molecules in the browser. Pick a molecule (or paste an XYZ file),
choose Hartree–Fock, LDA or PBE and a basis set, and see the electron density as a volumetric cloud
with molecular orbitals as isosurfaces. Drag an atom and the SCF re-runs as you move it; press
**Relax** and the geometry optimizes in front of you.

![Water HOMO, benzene HOMO and ethylene LUMO with the density cloud](docs/molecules.png)

## What it computes

| | |
|---|---|
| Methods | RHF, UHF (open shells), restricted Kohn–Sham DFT with LDA (Slater + VWN5) and PBE |
| Basis sets | STO-3G, 6-31G, 6-31G\* (cartesian d), elements H–Ne |
| Properties | total energy, orbital energies, HOMO–LUMO gap, dipole moment, analytic nuclear gradients |
| Geometry | BFGS optimization with analytic gradients, warm-started SCF between steps |
| Presets | H₂, H₂O, NH₃, CH₄, HF, N₂, CO₂, CH₂O, C₂H₄, C₆H₆ (STO-3G), OH radical, triplet O₂ |

Everything is checked against [PySCF](https://pyscf.org) (see `scripts/gen_reference.py`): energies
agree to 1e-8 Ha for Hartree–Fock and 2e-5 Ha for DFT, gradients to 1e-6 and 1e-4 Ha/bohr.

## How it works

```
crates/qchem        Rust library: integrals, SCF, DFT grids and functionals, gradients, BFGS
crates/qchem-wasm   wasm-bindgen bindings (Session, Optimizer)
src/engine          Web Worker hosting the WASM engine, and the client the UI talks to
src/render          WebGPU: compute shader for the density/orbital volume, ray marcher, atoms
src/ui, src/main.ts Controls, orbital level diagram, readouts, dragging
```

- **Quantum chemistry runs in Rust, compiled to WebAssembly**, in f64. Integrals use the
  McMurchie–Davidson scheme with Schwarz screening; the SCF starts from a superposition of atomic
  densities and uses DIIS. DFT uses Becke-partitioned Mura–Knowles/Lebedev grids. The library has no
  browser dependencies, so it can back other simulators or run natively
  (`cargo run --release --example water`).
- **The GPU does the volumes.** WGSL has no f64, but rendering does not need it: after each SCF the
  worker sends the basis and orbital coefficients, and a compute shader evaluates the density and
  the selected orbital on a ~112³ grid into an `rgba16float` 3D texture. A full-screen pass
  ray-marches it (emission/absorption for the cloud, sign-change detection for the ± isosurfaces)
  and stops at the depth of the ball-and-stick atoms.
- **SCF playback** records the orbitals of every iteration, so the cloud can be replayed from the
  atomic guess to convergence without recomputing.

Timings for one SCF on a single desktop core (measured in Node with the release WASM build):

| | Hartree–Fock | PBE |
|---|---|---|
| water, 6-31G\* | 15 ms | 0.3 s |
| ethylene, 6-31G\* | 0.13 s | 1.3 s |
| benzene, STO-3G | 0.6 s | 2.5 s |

DFT time is dominated by the integration grid, which is the next thing to move to the GPU.

## Running it

Needs Rust with the `wasm32-unknown-unknown` target, [wasm-pack](https://rustwasm.github.io/wasm-pack/)
and Node 20+.

```sh
npm install
npm run build:wasm   # compiles crates/qchem-wasm into crates/qchem-wasm/pkg
npm run dev          # http://localhost:5173
```

The 3D view needs WebGPU (recent Chrome, Edge or Safari). Without it the calculations still run and
the readouts work. `?res=64` lowers the volume resolution for slow GPUs.

## Tests

```sh
cargo test --release   # Rust: integrals, all methods vs PySCF, finite-difference gradients, invariances
npm test               # Vitest: the WASM build vs PySCF, worker protocol, presets, orbital evaluation
npm run test:e2e       # Playwright: the app in headless Chromium (WebGPU on SwiftShader)
```

The browser tests open the app with `?headless`, which renders offscreen because headless Chromium
cannot present a WebGPU canvas. They read the GPU volume back and compare it point by point with
the CPU evaluation, and exercise presets, method switching, dragging, relaxing and SCF playback.

To regenerate the PySCF reference data: `pip install pyscf && python scripts/gen_reference.py`.

## Not yet

Open-shell DFT, molecular dynamics, GPU two-electron integrals for larger molecules, density
differences and electrostatic potential maps.
