import { BOHR_PER_ANGSTROM } from './elements';
import type { Atom, Vec3 } from './molecule';

/**
 * A crystal is a lattice: one cell of atoms repeated forever in every direction. The quantum
 * chemistry here works on finite molecules (no periodic boundary conditions), so a crystal preset
 * cuts a cluster out of the lattice and computes that. Covalent crystals get hydrogen atoms on the
 * bonds that were cut, the usual way to make a cluster behave like the inside of the crystal
 * instead of a piece with broken bonds. The rest of the lattice is drawn around it, faded, and
 * is not part of the calculation.
 */
export interface CrystalSpec {
  /** Lattice vectors of the conventional cell (Angstrom). */
  cell: [Vec3, Vec3, Vec3];
  /** Atoms of one cell in fractional coordinates. */
  sites: { symbol: string; frac: Vec3 }[];
  /** Centre of the computed cluster, fractional. */
  center: Vec3;
  /** Lattice sites closer than this to the centre are computed (Angstrom). */
  radius: number;
  /** Nearest-neighbour bond length (Angstrom); bonds are drawn and cut up to 1.15 times this. */
  bond: number;
  /** Put hydrogen on cut bonds at this C-H, Si-H... distance (Angstrom). Not for ionic crystals. */
  capBond?: number;
  /** Cells repeated along each lattice vector when the preset opens. */
  repeats: number;
  /** Most cells along each lattice vector the slider offers. */
  maxRepeats: number;
  /** Lattice vectors that stay at one cell whatever `repeats` is (graphite's stacking axis). */
  fixedAxes?: number[];
}

export interface CrystalScene {
  /** The computed cluster, hydrogen caps last (bohr). */
  cluster: Atom[];
  /** How many of `cluster` are hydrogen caps. */
  caps: number;
  /** The lattice around it, drawn but not computed (bohr). */
  ghosts: Atom[];
  /**
   * Bonds from a ghost (index into `ghosts`) to another ghost or to a lattice atom of the
   * cluster (index into `cluster`), so they follow the cluster atoms when those move.
   */
  ghostBonds: { ghost: number; other: number; otherIsCluster: boolean }[];
  /** Edges of every repeated cell (bohr). */
  cellEdges: [Vec3, Vec3][];
}

const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a: Vec3, b: Vec3): Vec3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const scale = (a: Vec3, s: number): Vec3 => [a[0] * s, a[1] * s, a[2] * s];
const dist = (a: Vec3, b: Vec3) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

function cartesian(cell: CrystalSpec['cell'], f: Vec3): Vec3 {
  return add(add(scale(cell[0], f[0]), scale(cell[1], f[1])), scale(cell[2], f[2]));
}

/** Lattice sites in cells i, j, k in [lo, hi] (Angstrom), each once. */
function latticeSites(spec: CrystalSpec, lo: number, hi: number): { symbol: string; frac: Vec3; xyz: Vec3 }[] {
  const out: { symbol: string; frac: Vec3; xyz: Vec3 }[] = [];
  for (let i = lo; i <= hi; i++)
    for (let j = lo; j <= hi; j++)
      for (let k = lo; k <= hi; k++)
        for (const s of spec.sites) {
          const frac: Vec3 = [s.frac[0] + i, s.frac[1] + j, s.frac[2] + k];
          out.push({ symbol: s.symbol, frac, xyz: cartesian(spec.cell, frac) });
        }
  return out;
}

/**
 * The computed cluster and the lattice drawn around it, with `repeats` cells along each lattice
 * vector, centred on the cluster. Everything is shifted so the cluster centre is the origin.
 */
export function buildCrystal(spec: CrystalSpec, repeats = spec.repeats): CrystalScene {
  const center = cartesian(spec.cell, spec.center);
  const cutoff = spec.bond * 1.15;
  // Enough cells around the centre to hold the cluster and its neighbours.
  const reach = Math.ceil((spec.radius + 2 * spec.bond) / Math.min(...spec.cell.map((v) => Math.hypot(...v)))) + 1;
  const near = latticeSites(spec, -reach - 1, reach + 1);

  // The cluster: sites within the radius, then drop atoms left hanging by fewer than two bonds
  // (a lone CH3 group or an isolated ion at a corner is not a piece of the crystal).
  let inside = near.filter((s) => dist(s.xyz, center) < spec.radius);
  for (;;) {
    const kept = inside.filter((s) => inside.filter((t) => t !== s && dist(s.xyz, t.xyz) < cutoff).length >= 2);
    if (kept.length === inside.length) break;
    inside = kept;
  }
  const caps: Atom[] = [];
  if (spec.capBond !== undefined) {
    for (const s of inside)
      for (const t of near)
        if (dist(s.xyz, t.xyz) < cutoff && dist(s.xyz, t.xyz) > 0.1 && !inside.includes(t)) {
          const dir = scale(sub(t.xyz, s.xyz), 1 / dist(s.xyz, t.xyz));
          caps.push({ symbol: 'H', xyz: add(s.xyz, scale(dir, spec.capBond)) });
        }
  }
  const toBohr = (p: Vec3): Vec3 => scale(sub(p, center), BOHR_PER_ANGSTROM);
  const cluster = [...inside.map((s) => ({ symbol: s.symbol, xyz: toBohr(s.xyz) })), ...caps.map((c) => ({ symbol: 'H', xyz: toBohr(c.xyz) }))];

  // Repeated cells centred on the cluster: cells n0 .. n0 + n - 1 along each vector.
  const n = [0, 1, 2].map((d) => (spec.fixedAxes?.includes(d) ? 1 : repeats));
  const n0 = spec.center.map((c, d) => Math.floor(c - n[d] / 2 + 0.5)) as Vec3;
  const eps = 1e-6;
  const block = latticeSites(spec, Math.min(...n0) - 1, Math.max(...n0) + repeats + 1).filter((s) =>
    s.frac.every((f, d) => f >= n0[d] - eps && f <= n0[d] + n[d] + eps),
  );
  const isCluster = (p: Vec3) => inside.some((s) => dist(s.xyz, p) < 1e-3);
  const ghostSites = block.filter((s) => !isCluster(s.xyz));
  const ghosts = ghostSites.map((s) => ({ symbol: s.symbol, xyz: toBohr(s.xyz) }));

  // Bonds from each ghost to its neighbours among the ghosts and the cluster's lattice atoms
  // (not its hydrogen caps, which stand in for exactly those neighbours).
  const ghostBonds: CrystalScene['ghostBonds'] = [];
  const bonded = (a: Vec3, b: Vec3) => {
    const d = dist(a, b);
    return d < cutoff && d > 0.1;
  };
  ghostSites.forEach((g, i) => {
    for (let j = i + 1; j < ghostSites.length; j++)
      if (bonded(g.xyz, ghostSites[j].xyz)) ghostBonds.push({ ghost: i, other: j, otherIsCluster: false });
    inside.forEach((c, k) => {
      if (bonded(g.xyz, c.xyz)) ghostBonds.push({ ghost: i, other: k, otherIsCluster: true });
    });
  });

  const cellEdges: [Vec3, Vec3][] = [];
  const corner = (i: number, j: number, k: number) => toBohr(cartesian(spec.cell, [n0[0] + i, n0[1] + j, n0[2] + k]));
  for (let i = 0; i <= n[0]; i++)
    for (let j = 0; j <= n[1]; j++)
      for (let k = 0; k <= n[2]; k++) {
        if (i < n[0]) cellEdges.push([corner(i, j, k), corner(i + 1, j, k)]);
        if (j < n[1]) cellEdges.push([corner(i, j, k), corner(i, j + 1, k)]);
        if (k < n[2]) cellEdges.push([corner(i, j, k), corner(i, j, k + 1)]);
      }
  return { cluster, caps: caps.length, ghosts, ghostBonds, cellEdges };
}

const DIAMOND_SITES = (symbol: string) =>
  (
    [
      [0, 0, 0],
      [0, 0.5, 0.5],
      [0.5, 0, 0.5],
      [0.5, 0.5, 0],
      [0.25, 0.25, 0.25],
      [0.25, 0.75, 0.75],
      [0.75, 0.25, 0.75],
      [0.75, 0.75, 0.25],
    ] as Vec3[]
  ).map((frac) => ({ symbol, frac }));

const ROCK_SALT_SITES = (cation: string, anion: string) =>
  (
    [
      [0, 0, 0],
      [0, 0.5, 0.5],
      [0.5, 0, 0.5],
      [0.5, 0.5, 0],
    ] as Vec3[]
  ).flatMap((f) => [
    { symbol: cation, frac: f },
    { symbol: anion, frac: [f[0] + 0.5, f[1], f[2]] as Vec3 },
  ]);

const cubic = (a: number): [Vec3, Vec3, Vec3] => [
  [a, 0, 0],
  [0, a, 0],
  [0, 0, a],
];

/** Experimental lattice constants at room temperature. */
export const CRYSTALS: Record<string, CrystalSpec> = {
  // The ten carbons nearest the centre of the cubic cell form one adamantane cage.
  diamond: {
    cell: cubic(3.567),
    sites: DIAMOND_SITES('C'),
    center: [0.5, 0.5, 0.5],
    radius: 1.9,
    bond: 1.5445,
    capBond: 1.09,
    repeats: 2,
    maxRepeats: 3,
  },
  silicon: {
    cell: cubic(5.431),
    sites: DIAMOND_SITES('Si'),
    center: [0.5, 0.5, 0.5],
    radius: 2.9,
    bond: 2.3517,
    capBond: 1.48,
    repeats: 2,
    maxRepeats: 3,
  },
  // Eight ions on the corners of one small cube, an eighth of the cubic cell.
  nacl: {
    cell: cubic(5.64),
    sites: ROCK_SALT_SITES('Na', 'Cl'),
    center: [0.25, 0.25, 0.25],
    radius: 2.5,
    bond: 2.82,
    repeats: 1,
    maxRepeats: 3,
  },
  mgo: {
    cell: cubic(4.212),
    sites: ROCK_SALT_SITES('Mg', 'O'),
    center: [0.25, 0.25, 0.25],
    radius: 1.9,
    bond: 2.106,
    repeats: 1,
    maxRepeats: 4,
  },
  // One graphene sheet of graphite (AB stacking); the cluster is a pyrene-shaped flake.
  graphite: {
    cell: [
      [2.464, 0, 0],
      [-1.232, 2.1339, 0],
      [0, 0, 6.711],
    ],
    sites: [
      { symbol: 'C', frac: [0, 0, 0.25] },
      { symbol: 'C', frac: [1 / 3, 2 / 3, 0.25] },
      { symbol: 'C', frac: [0, 0, 0.75] },
      { symbol: 'C', frac: [2 / 3, 1 / 3, 0.75] },
    ],
    center: [1 / 6, 1 / 3, 0.25],
    radius: 3.8,
    bond: 1.4226,
    capBond: 1.08,
    repeats: 4,
    maxRepeats: 6,
    fixedAxes: [2],
  },
};
