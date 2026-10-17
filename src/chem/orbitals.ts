import type { ComputeOutput, GpuBasis } from '../engine/types';

/**
 * CPU version of the basis-function evaluation done in `volume.wgsl`. It is the reference the
 * shader is checked against, and the GpuBasis layout is defined by `gpu_basis` in the WASM crate.
 */
export function evalAo(basis: GpuBasis, x: number, y: number, z: number, out: Float64Array): void {
  const { shells, prims, componentNorms } = basis;
  const nshell = shells.length / 8;
  for (let s = 0; s < nshell; s++) {
    const o = 8 * s;
    const dx = x - shells[o];
    const dy = y - shells[o + 1];
    const dz = z - shells[o + 2];
    const l = shells[o + 3];
    const nprim = shells[o + 4];
    const p0 = shells[o + 5];
    const f = shells[o + 6];
    const r2 = dx * dx + dy * dy + dz * dz;
    let radial = 0;
    for (let p = 0; p < nprim; p++) radial += prims[2 * (p0 + p) + 1] * Math.exp(-prims[2 * (p0 + p)] * r2);
    const n = componentNorms[l];
    if (l === 0) out[f] = radial * n[0];
    else if (l === 1) {
      out[f] = radial * dx * n[0];
      out[f + 1] = radial * dy * n[1];
      out[f + 2] = radial * dz * n[2];
    } else {
      out[f] = radial * dx * dx * n[0];
      out[f + 1] = radial * dx * dy * n[1];
      out[f + 2] = radial * dx * dz * n[2];
      out[f + 3] = radial * dy * dy * n[3];
      out[f + 4] = radial * dy * dz * n[4];
      out[f + 5] = radial * dz * dz * n[5];
    }
  }
}

/** Value of molecular orbital `mo` (columns of an nbf x nmo row-major matrix) from AO values. */
export function moValue(ao: Float64Array, c: ArrayLike<number>, nmo: number, mo: number): number {
  let v = 0;
  for (let i = 0; i < ao.length; i++) v += ao[i] * c[i * nmo + mo];
  return v;
}

/** Electron density at a point from the occupied alpha and beta orbitals. */
export function densityAt(res: ComputeOutput, x: number, y: number, z: number, ao = new Float64Array(res.nbf)): number {
  evalAo(res.basis, x, y, z, ao);
  let rho = 0;
  for (let k = 0; k < res.nalpha; k++) rho += moValue(ao, res.ca, res.nmo, k) ** 2;
  const cb = res.unrestricted ? res.cb : res.ca;
  for (let k = 0; k < res.nbeta; k++) rho += moValue(ao, cb, res.nmo, k) ** 2;
  return rho;
}

export interface OrbitalLevel {
  index: number;
  energy: number;
  occupied: boolean;
  label: string;
}

/** Orbital levels with HOMO/LUMO labels (alpha spin for open shells). */
export function orbitalLevels(res: ComputeOutput, spin: 'alpha' | 'beta' = 'alpha'): OrbitalLevel[] {
  const eps = spin === 'beta' && res.unrestricted ? res.epsB : res.epsA;
  const nocc = spin === 'beta' ? res.nbeta : res.nalpha;
  return eps.map((energy, index) => {
    let label = '';
    if (index === nocc - 1) label = 'HOMO';
    else if (index === nocc) label = 'LUMO';
    else if (index < nocc) label = `HOMO−${nocc - 1 - index}`;
    else label = `LUMO+${index - nocc}`;
    return { index, energy, occupied: index < nocc, label };
  });
}

/** HOMO and LUMO energies (Hartree). For open shells, the highest occupied and lowest empty over both spins. */
export function frontier(res: ComputeOutput): { homo: number; lumo: number | null } {
  let homo = res.epsA[res.nalpha - 1] ?? -Infinity;
  let lumo = res.epsA[res.nalpha] ?? Infinity;
  if (res.unrestricted) {
    if (res.nbeta > 0) homo = Math.max(homo, res.epsB[res.nbeta - 1]);
    lumo = Math.min(lumo, res.epsB[res.nbeta] ?? Infinity);
  }
  return { homo, lumo: Number.isFinite(lumo) ? lumo : null };
}
