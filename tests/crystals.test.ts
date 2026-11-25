import { describe, expect, it } from 'vitest';
import { CRYSTALS, buildCrystal } from '../src/chem/crystals';
import { BOHR_PER_ANGSTROM } from '../src/chem/elements';
import type { Atom } from '../src/chem/molecule';

const dist = (a: Atom, b: Atom) => Math.hypot(...[0, 1, 2].map((k) => a.xyz[k] - b.xyz[k])) / BOHR_PER_ANGSTROM;
const formula = (atoms: Atom[]) => {
  const n: Record<string, number> = {};
  for (const a of atoms) n[a.symbol] = (n[a.symbol] ?? 0) + 1;
  return n;
};

describe('crystal pieces', () => {
  it('cut the expected clusters', () => {
    expect(formula(buildCrystal(CRYSTALS.diamond).cluster)).toEqual({ C: 10, H: 16 });
    expect(formula(buildCrystal(CRYSTALS.silicon).cluster)).toEqual({ Si: 10, H: 16 });
    expect(formula(buildCrystal(CRYSTALS.graphite).cluster)).toEqual({ C: 16, H: 10 });
    expect(formula(buildCrystal(CRYSTALS.nacl).cluster)).toEqual({ Na: 4, Cl: 4 });
    expect(formula(buildCrystal(CRYSTALS.mgo).cluster)).toEqual({ Mg: 4, O: 4 });
  });

  it('give every diamond carbon four neighbours: carbons or hydrogen caps', () => {
    const { cluster } = buildCrystal(CRYSTALS.diamond);
    for (const c of cluster.filter((a) => a.symbol === 'C')) {
      const cc = cluster.filter((b) => b.symbol === 'C' && b !== c && Math.abs(dist(b, c) - 1.5445) < 1e-3);
      const ch = cluster.filter((b) => b.symbol === 'H' && Math.abs(dist(b, c) - 1.09) < 1e-3);
      expect(cc.length + ch.length).toBe(4);
    }
  });

  it('draw the rest of the lattice without doubling the computed atoms', () => {
    for (const spec of Object.values(CRYSTALS)) {
      for (let n = 1; n <= spec.maxRepeats; n++) {
        const s = buildCrystal(spec, [n, n, n]);
        expect(s.ghosts.length).toBeGreaterThan(0);
        expect(s.ghosts.every((g) => g.symbol !== 'H')).toBe(true);
        const lattice = s.cluster.slice(0, s.cluster.length - s.caps);
        for (const g of s.ghosts) for (const c of lattice) expect(dist(g, c)).toBeGreaterThan(0.5);
        for (const b of s.ghostBonds) expect(b.ghost).toBeLessThan(s.ghosts.length);
      }
    }
  });

  it('use the smallest repeating cell', () => {
    // Face-centred cubic: two atoms per cell, a quarter of the cube of side a.
    for (const id of ['diamond', 'silicon', 'nacl', 'mgo']) {
      const [a, b, c] = CRYSTALS[id].cell;
      const volume = Math.abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]));
      const side = Math.SQRT2 * Math.hypot(...a);
      expect(CRYSTALS[id].sites, id).toHaveLength(2);
      expect(volume / side ** 3, id).toBeCloseTo(1 / 4, 10);
    }
  });

  it('grow the drawn block with Nx, Ny and Nz', () => {
    const edges = ([x, y, z]: number[]) => x * (y + 1) * (z + 1) + y * (x + 1) * (z + 1) + z * (x + 1) * (y + 1);
    for (const n of [[1, 1, 1], [2, 3, 1], [3, 3, 3], [1, 4, 2]] as [number, number, number][]) {
      expect(buildCrystal(CRYSTALS.nacl, n).cellEdges.length).toBe(edges(n));
      expect(buildCrystal(CRYSTALS.graphite, n).cellEdges.length).toBe(edges(n));
    }
  });
});
