import { BOHR_PER_ANGSTROM, element } from '../chem/elements';
import { bonds, type Atom } from '../chem/molecule';
import type { BondInstance, SphereInstance } from '../render/renderer';

export const BOND_RADIUS = 0.17;

/** Ball-and-stick radius in bohr. */
export function atomRadius(symbol: string): number {
  return 0.35 + 0.3 * element(symbol).covalentRadius * BOHR_PER_ANGSTROM;
}

/** Display colours tuned for a dark background. */
export function atomColor(symbol: string): [number, number, number] {
  return element(symbol).color;
}

export function buildScene(atoms: Atom[], highlight: number | null): { spheres: SphereInstance[]; bonds: BondInstance[] } {
  const spheres = atoms.map((a, i) => ({
    center: a.xyz,
    radius: atomRadius(a.symbol),
    color: atomColor(a.symbol),
    highlight: i === highlight,
  }));
  const sticks = bonds(atoms).map(([i, j]) => ({
    a: atoms[i].xyz,
    b: atoms[j].xyz,
    radius: BOND_RADIUS,
    colorA: atomColor(atoms[i].symbol),
    colorB: atomColor(atoms[j].symbol),
  }));
  return { spheres, bonds: sticks };
}

/** Bounding sphere of the nuclei (bohr). */
export function boundingSphere(atoms: Atom[]): { center: [number, number, number]; radius: number } {
  const c = [0, 1, 2].map((k) => atoms.reduce((s, a) => s + a.xyz[k], 0) / Math.max(1, atoms.length)) as [number, number, number];
  const radius = atoms.reduce((r, a) => Math.max(r, Math.hypot(a.xyz[0] - c[0], a.xyz[1] - c[1], a.xyz[2] - c[2])), 0);
  return { center: c, radius };
}
