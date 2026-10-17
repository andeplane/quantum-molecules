import { DEBYE_PER_AU, EV_PER_HARTREE } from '../chem/elements';
import { frontier } from '../chem/orbitals';
import type { ComputeOutput } from '../engine/types';

export interface Readout {
  label: string;
  value: string;
  title?: string;
}

const ev = (h: number) => `${(h * EV_PER_HARTREE).toFixed(2)} eV`;

export function readouts(res: ComputeOutput): Readout[] {
  const { homo, lumo } = frontier(res);
  const dipole = Math.hypot(...res.dipole) * DEBYE_PER_AU;
  const rows: Readout[] = [
    { label: 'Total energy', value: `${res.energy.toFixed(6)} Ha`, title: `${(res.energy * EV_PER_HARTREE).toFixed(3)} eV` },
    { label: 'HOMO', value: ev(homo) },
    { label: 'LUMO', value: lumo === null ? 'none' : ev(lumo) },
    { label: 'HOMO–LUMO gap', value: lumo === null ? 'n/a' : ev(lumo - homo) },
    { label: 'Dipole moment', value: `${dipole.toFixed(3)} D` },
    {
      label: 'SCF',
      value: `${res.converged ? 'converged' : 'NOT converged'} in ${res.iterations} iterations, ${Math.round(res.millis)} ms`,
    },
    { label: 'Basis functions', value: String(res.nbf) },
  ];
  if (res.maxForce !== undefined) rows.push({ label: 'Max force', value: `${res.maxForce.toExponential(2)} Ha/bohr` });
  return rows;
}
