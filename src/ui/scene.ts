import type { CrystalScene } from '../chem/crystals';
import { BOHR_PER_ANGSTROM, element } from '../chem/elements';
import { bonds, type Atom } from '../chem/molecule';
import type { BondInstance, SphereInstance } from '../render/renderer';

export const BOND_RADIUS = 0.17;

/** Ball-and-stick radius in bohr. */
export function atomRadius(symbol: string): number {
  // Capped so the big alkali and alkaline-earth atoms (Na, Mg) do not swamp their neighbours.
  return 0.35 + 0.3 * Math.min(element(symbol).covalentRadius, 1.15) * BOHR_PER_ANGSTROM;
}

/** Display colours tuned for a dark background. */
export function atomColor(symbol: string): [number, number, number] {
  return element(symbol).color;
}

/**
 * Ball-and-stick instances. `maxBond` (Angstrom) caps the bond length, for crystal pieces where
 * the covalent-radius rule would join second neighbours such as Na and Na in salt.
 */
export function buildScene(
  atoms: Atom[],
  highlight: number | null,
  maxBond?: number,
): { spheres: SphereInstance[]; bonds: BondInstance[] } {
  const spheres = atoms.map((a, i) => ({
    center: a.xyz,
    radius: atomRadius(a.symbol),
    color: atomColor(a.symbol),
    highlight: i === highlight,
  }));
  const within = ([i, j]: [number, number]) =>
    maxBond === undefined || Math.hypot(...[0, 1, 2].map((k) => atoms[i].xyz[k] - atoms[j].xyz[k])) < maxBond * BOHR_PER_ANGSTROM;
  const sticks = bonds(atoms).filter(within).map(([i, j]) => ({
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

/** Unit-cell edges: a muted gold that reads as a guide, not as part of the molecule. */
export const CELL_COLOR: [number, number, number] = [0.79, 0.65, 0.35];

/**
 * The crystal around a computed cluster: faded lattice atoms and bonds, and the unit-cell edges.
 * Bonds to the cluster follow `atoms`, the cluster as it is now (its lattice atoms come first).
 */
export function crystalExtras(scene: CrystalScene, atoms: Atom[]): { spheres: SphereInstance[]; bonds: BondInstance[] } {
  const spheres: SphereInstance[] = scene.ghosts.map((a) => ({
    center: a.xyz,
    radius: atomRadius(a.symbol) * 0.6,
    color: atomColor(a.symbol),
    highlight: false,
    ghost: true,
  }));
  const sticks: BondInstance[] = scene.ghostBonds.map(({ ghost, other, otherIsCluster }) => {
    const a = scene.ghosts[ghost];
    const b = otherIsCluster ? atoms[other] : scene.ghosts[other];
    return { a: a.xyz, b: b.xyz, radius: BOND_RADIUS * 0.45, colorA: atomColor(a.symbol), colorB: atomColor(b.symbol), ghost: true };
  });
  for (const [a, b] of scene.cellEdges) sticks.push({ a, b, radius: 0.045, colorA: CELL_COLOR, colorB: CELL_COLOR });
  return { spheres, bonds: sticks };
}
