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
        const s = buildCrystal(spec, n);
        expect(s.ghosts.length).toBeGreaterThan(0);
        expect(s.ghosts.every((g) => g.symbol !== 'H')).toBe(true);
        const lattice = s.cluster.slice(0, s.cluster.length - s.caps);
        for (const g of s.ghosts) for (const c of lattice) expect(dist(g, c)).toBeGreaterThan(0.5);
        for (const b of s.ghostBonds) expect(b.ghost).toBeLessThan(s.ghosts.length);
      }
    }
  });

  it('grow the drawn block with the number of cells', () => {
    const cells = (n: number) => buildCrystal(CRYSTALS.nacl, n).cellEdges.length;
    // A block of n^3 cubes has 3 n (n + 1)^2 edges.
    expect([1, 2, 3].map(cells)).toEqual([12, 54, 144]);
    // Graphite repeats in the sheet only.
    expect(buildCrystal(CRYSTALS.graphite, 2).cellEdges.length).toBe(2 * (2 * 3 * 2) + 3 * 3);
  });
});
