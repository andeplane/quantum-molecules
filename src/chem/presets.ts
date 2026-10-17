import type { BasisName } from '../engine/types';

export interface Preset {
  id: string;
  name: string;
  /** XYZ body in Angstrom (symbol x y z per line). */
  xyz: string;
  charge: number;
  multiplicity: number;
  /** Basis sets that stay interactive for this molecule. */
  bases: BasisName[];
}

const ALL: BasisName[] = ['sto-3g', '6-31g', '6-31g*'];

function benzene(): string {
  const lines: string[] = [];
  for (let k = 0; k < 6; k++) {
    const t = (k * Math.PI) / 3;
    lines.push(`C ${(1.397 * Math.cos(t)).toFixed(4)} ${(1.397 * Math.sin(t)).toFixed(4)} 0`);
  }
  for (let k = 0; k < 6; k++) {
    const t = (k * Math.PI) / 3;
    lines.push(`H ${(2.481 * Math.cos(t)).toFixed(4)} ${(2.481 * Math.sin(t)).toFixed(4)} 0`);
  }
  return lines.join('\n');
}

/** Experimental equilibrium geometries (rounded). */
export const PRESETS: Preset[] = [
  { id: 'h2', name: 'Hydrogen H₂', xyz: 'H 0 0 -0.3705\nH 0 0 0.3705', charge: 0, multiplicity: 1, bases: ALL },
  {
    id: 'h2o',
    name: 'Water H₂O',
    xyz: 'O 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'nh3',
    name: 'Ammonia NH₃',
    xyz: 'N 0 0 0.1162\nH 0 0.9377 -0.2711\nH 0.8121 -0.4689 -0.2711\nH -0.8121 -0.4689 -0.2711',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'ch4',
    name: 'Methane CH₄',
    xyz: 'C 0 0 0\nH 0.6291 0.6291 0.6291\nH -0.6291 -0.6291 0.6291\nH -0.6291 0.6291 -0.6291\nH 0.6291 -0.6291 -0.6291',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  { id: 'hf', name: 'Hydrogen fluoride HF', xyz: 'F 0 0 0.0917\nH 0 0 -0.8253', charge: 0, multiplicity: 1, bases: ALL },
  { id: 'n2', name: 'Nitrogen N₂', xyz: 'N 0 0 -0.5488\nN 0 0 0.5488', charge: 0, multiplicity: 1, bases: ALL },
  { id: 'co2', name: 'Carbon dioxide CO₂', xyz: 'C 0 0 0\nO 0 0 1.1621\nO 0 0 -1.1621', charge: 0, multiplicity: 1, bases: ALL },
  {
    id: 'ch2o',
    name: 'Formaldehyde CH₂O',
    xyz: 'C 0 0 -0.5285\nO 0 0 0.6745\nH 0 0.9363 -1.1079\nH 0 -0.9363 -1.1079',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'c2h4',
    name: 'Ethylene C₂H₄',
    xyz: 'C 0 0 0.6695\nC 0 0 -0.6695\nH 0 0.9288 1.2321\nH 0 -0.9288 1.2321\nH 0 0.9288 -1.2321\nH 0 -0.9288 -1.2321',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  { id: 'c6h6', name: 'Benzene C₆H₆', xyz: benzene(), charge: 0, multiplicity: 1, bases: ['sto-3g'] },
  { id: 'oh', name: 'Hydroxyl radical OH', xyz: 'O 0 0 0\nH 0 0 0.9697', charge: 0, multiplicity: 2, bases: ALL },
  { id: 'o2', name: 'Oxygen O₂ (triplet)', xyz: 'O 0 0 -0.6038\nO 0 0 0.6038', charge: 0, multiplicity: 3, bases: ALL },
];

export function preset(id: string): Preset {
  const p = PRESETS.find((p) => p.id === id);
  if (!p) throw new Error(`Unknown preset ${id}`);
  return p;
}
