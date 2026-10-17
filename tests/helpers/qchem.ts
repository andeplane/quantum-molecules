import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import * as qchem from '@qchem/qchem_wasm.js';
import type { ComputeInput, ComputeOutput } from '../../src/engine/types';

const wasm = fileURLToPath(new URL('../../crates/qchem-wasm/pkg/qchem_wasm_bg.wasm', import.meta.url));
qchem.initSync({ module: readFileSync(wasm) });

export { qchem };

export interface RefCase {
  name: string;
  basis: ComputeInput['basis'];
  method: 'rhf' | 'uhf' | 'rks';
  xc: string | null;
  charge: number;
  multiplicity: number;
  atoms: { symbol: string; xyz: [number, number, number] }[];
  nbf: number;
  energy: number;
  dipole: [number, number, number];
  moEnergies?: number[];
  gradient: [number, number, number][];
}

const fixtures = fileURLToPath(new URL('../../crates/qchem/tests/fixtures/pyscf-reference.json', import.meta.url));
export const REFERENCE: RefCase[] = JSON.parse(readFileSync(fixtures, 'utf8'));

export function refCase(name: string, basis: string, xc: string | null = null): RefCase {
  const c = REFERENCE.find((c) => c.name === name && c.basis === basis && c.xc === xc);
  if (!c) throw new Error(`no reference case ${name}/${basis}/${xc}`);
  return c;
}

export function inputFor(c: RefCase, extra: Partial<ComputeInput> = {}): ComputeInput {
  const method = c.xc === null ? 'hf' : c.xc.startsWith('pbe') ? 'pbe' : 'lda';
  return { atoms: c.atoms, charge: c.charge, multiplicity: c.multiplicity, basis: c.basis, method, ...extra };
}

export function compute(input: ComputeInput, session = new qchem.Session()): ComputeOutput {
  return session.compute(input) as ComputeOutput;
}
