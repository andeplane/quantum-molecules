import { element } from '../chem/elements';
import type { Atom } from '../chem/molecule';

const NAMES: Record<string, string> = {
  H: 'hydrogen',
  He: 'helium',
  Li: 'lithium',
  Be: 'beryllium',
  B: 'boron',
  C: 'carbon',
  N: 'nitrogen',
  O: 'oxygen',
  F: 'fluorine',
  Ne: 'neon',
};

/** Orbital lobe colours; keep in sync with `raymarch.wgsl`. */
export const LOBE_POSITIVE: [number, number, number] = [0.22, 0.52, 1.0];
export const LOBE_NEGATIVE: [number, number, number] = [0.95, 0.42, 0.18];
/** Stops of the density cloud ramp; keep in sync with `cloudColor` in `raymarch.wgsl`. */
export const CLOUD_RAMP: [number, number, number][] = [
  [0.3, 0.28, 0.4],
  [0.66, 0.64, 0.78],
  [1.0, 0.98, 0.95],
];

export const css = ([r, g, b]: [number, number, number]) => `rgb(${Math.round(r * 255)}, ${Math.round(g * 255)}, ${Math.round(b * 255)})`;

export interface LegendElement {
  symbol: string;
  name: string;
  color: string;
}

/** The elements in a molecule, heaviest first, with their display colours. */
export function legendElements(atoms: Atom[]): LegendElement[] {
  const symbols = [...new Set(atoms.map((a) => element(a.symbol).symbol))];
  symbols.sort((a, b) => element(b).z - element(a).z);
  return symbols.map((symbol) => ({ symbol, name: NAMES[symbol] ?? symbol, color: css(element(symbol).color) }));
}

/** A plain-language description of an orbital label such as HOMO or LUMO+1. */
export function describeOrbital(label: string): string {
  if (label === 'HOMO') return 'the highest-energy orbital that holds electrons';
  if (label === 'LUMO') return 'the lowest-energy empty orbital';
  if (label.startsWith('HOMO')) return 'an occupied orbital below the HOMO';
  return 'an empty orbital above the LUMO';
}
