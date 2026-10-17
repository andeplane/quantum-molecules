export const BOHR_PER_ANGSTROM = 1 / 0.52917721092;
export const EV_PER_HARTREE = 27.211386245988;
export const DEBYE_PER_AU = 2.541746473;

export interface ElementInfo {
  symbol: string;
  z: number;
  /** Atomic mass in amu. */
  mass: number;
  /** Covalent radius in Angstrom. */
  covalentRadius: number;
  /** Display colour (CPK-like), as linear RGB 0..1. */
  color: [number, number, number];
}

const RAW: [string, number, number, [number, number, number]][] = [
  ['H', 1.008, 0.31, [0.95, 0.95, 0.95]],
  ['He', 4.0026, 0.28, [0.85, 1.0, 1.0]],
  ['Li', 6.94, 1.28, [0.8, 0.5, 1.0]],
  ['Be', 9.0122, 0.96, [0.76, 1.0, 0.0]],
  ['B', 10.81, 0.84, [1.0, 0.71, 0.71]],
  ['C', 12.011, 0.76, [0.38, 0.38, 0.38]],
  ['N', 14.007, 0.71, [0.19, 0.31, 0.97]],
  ['O', 15.999, 0.66, [1.0, 0.13, 0.13]],
  ['F', 18.998, 0.57, [0.56, 0.88, 0.31]],
  ['Ne', 20.18, 0.58, [0.7, 0.89, 0.96]],
];

export const ELEMENTS: ElementInfo[] = RAW.map(([symbol, mass, covalentRadius, color], i) => ({
  symbol,
  z: i + 1,
  mass,
  covalentRadius,
  color,
}));

const BY_SYMBOL = new Map(ELEMENTS.map((e) => [e.symbol.toLowerCase(), e]));

export function element(symbol: string): ElementInfo {
  const e = BY_SYMBOL.get(symbol.toLowerCase());
  if (!e) throw new Error(`Unsupported element: ${symbol} (only H to Ne are available)`);
  return e;
}
