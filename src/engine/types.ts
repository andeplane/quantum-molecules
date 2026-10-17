/** Mirrors `ComputeInput` / `ComputeOutput` in crates/qchem-wasm/src/lib.rs. */

export type BasisName = 'sto-3g' | '6-31g' | '6-31g*';
export type MethodName = 'hf' | 'lda' | 'pbe';
export type GridName = 'coarse' | 'medium' | 'fine' | 'ultrafine';

export interface AtomInput {
  symbol: string;
  /** Bohr. */
  xyz: [number, number, number];
}

export interface ComputeInput {
  atoms: AtomInput[];
  charge: number;
  multiplicity: number;
  basis: BasisName;
  method: MethodName;
  grid?: GridName;
  gradient?: boolean;
  history?: boolean;
  warmStart?: boolean;
}

export interface GpuBasis {
  /** Per shell: cx, cy, cz, l, nprim, primStart, fnOffset, 0. */
  shells: number[];
  /** Per primitive: exponent, coefficient. */
  prims: number[];
  componentNorms: number[][];
}

export interface IterationOutput {
  energy: number;
  error: number;
  ca: number[];
  cb: number[];
}

export interface ComputeOutput {
  energy: number;
  converged: boolean;
  iterations: number;
  unrestricted: boolean;
  nbf: number;
  nmo: number;
  nalpha: number;
  nbeta: number;
  epsA: number[];
  epsB: number[];
  ca: number[];
  cb: number[];
  dipole: [number, number, number];
  nuclearRepulsion: number;
  exc: number;
  gradient?: number[];
  maxForce?: number;
  history: IterationOutput[];
  basis: GpuBasis;
  coordinates: number[];
  millis: number;
}

export type WorkerRequest =
  | { id: number; type: 'compute'; input: ComputeInput }
  | { id: number; type: 'relax'; input: ComputeInput; maxSteps?: number }
  | { id: number; type: 'cancel' }
  | { id: number; type: 'reset' };

export type WorkerResponse =
  | { id: number; type: 'ready' }
  | { id: number; type: 'result'; output: ComputeOutput }
  | { id: number; type: 'relax-step'; output: ComputeOutput; step: number }
  | { id: number; type: 'relax-done'; converged: boolean; steps: number }
  | { id: number; type: 'error'; message: string };
