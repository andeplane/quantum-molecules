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
    {
      label: 'Total energy',
      value: `${res.energy.toFixed(6)} Ha`,
      title: `${(res.energy * EV_PER_HARTREE).toFixed(3)} eV. The energy of the electrons plus the repulsion between the nuclei; lower is more stable (1 Ha = 27.211 eV).`,
    },
    { label: 'HOMO', value: ev(homo), title: 'Energy of the highest occupied orbital. Its negative is roughly the energy needed to remove an electron.' },
    { label: 'LUMO', value: lumo === null ? 'none' : ev(lumo), title: 'Energy of the lowest empty orbital, where an added electron would go.' },
    { label: 'HOMO–LUMO gap', value: lumo === null ? 'n/a' : ev(lumo - homo), title: 'Roughly how much energy it takes to excite the molecule.' },
    { label: 'Dipole moment', value: `${dipole.toFixed(3)} D`, title: 'How unevenly the charge is spread, in debye. Water is about 1.85 D in experiment.' },
    {
      label: 'SCF',
      value: `${res.converged ? 'converged' : 'NOT converged'} in ${res.iterations} iterations, ${Math.round(res.millis)} ms`,
      title: 'The self-consistent field loop: guess the orbitals, compute the field they make, solve again, until nothing changes.',
    },
    { label: 'Basis functions', value: String(res.nbf), title: 'How many atom-centred functions the orbitals are built from. Bigger basis sets are more accurate and slower.' },
  ];
  if (res.maxForce !== undefined) rows.push({ label: 'Max force', value: `${res.maxForce.toExponential(2)} Ha/bohr` });
  return rows;
}
