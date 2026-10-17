import { describe, expect, it } from 'vitest';
import { BOHR_PER_ANGSTROM, ELEMENTS, element } from '../src/chem/elements';
import { bonds, centered, electronCount, parseXyz, toXyz } from '../src/chem/molecule';
import { PRESETS, preset } from '../src/chem/presets';
import { densityAt, evalAo, frontier, moValue, orbitalLevels } from '../src/chem/orbitals';
import type { ComputeOutput } from '../src/engine/types';
import { compute, inputFor, qchem, refCase } from './helpers/qchem';

describe('elements', () => {
  it('covers H to Ne in order', () => {
    expect(ELEMENTS.map((e) => e.symbol)).toEqual(['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']);
    expect(element('o').z).toBe(8);
    expect(() => element('Na')).toThrow(/Na/);
  });
});

describe('XYZ parsing', () => {
  it('reads a standard XYZ file in Angstrom', () => {
    const atoms = parseXyz('3\nwater\nO 0 0 0.1173\nh 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n');
    expect(atoms.map((a) => a.symbol)).toEqual(['O', 'H', 'H']);
    expect(atoms[1].xyz[1]).toBeCloseTo(0.7572 * BOHR_PER_ANGSTROM, 12);
  });

  it('skips lines that are not atoms and rejects unknown elements', () => {
    expect(parseXyz('# comment\n\nC 0 0 0 extra\nnot an atom line here')).toHaveLength(1);
    expect(() => parseXyz('U 0 0 0')).toThrow(/U/);
  });

  it('round-trips through toXyz', () => {
    const atoms = parseXyz(preset('ch2o').xyz);
    const back = parseXyz(toXyz(atoms, 'formaldehyde'));
    expect(toXyz(atoms).split('\n')[0]).toBe('4');
    back.forEach((a, i) => a.xyz.forEach((c, k) => expect(c).toBeCloseTo(atoms[i].xyz[k], 5)));
  });
});

describe('molecule helpers', () => {
  it('finds bonds from covalent radii', () => {
    expect(bonds(parseXyz(preset('h2o').xyz))).toEqual([
      [0, 1],
      [0, 2],
    ]);
    expect(bonds(parseXyz(preset('c6h6').xyz))).toHaveLength(12);
    expect(bonds(parseXyz('H 0 0 0\nH 0 0 3'))).toEqual([]);
  });

  it('centres molecules', () => {
    const c = centered(parseXyz(preset('h2o').xyz));
    for (let k = 0; k < 3; k++) expect(c.reduce((s, a) => s + a.xyz[k], 0)).toBeCloseTo(0, 12);
  });

  it('counts electrons', () => {
    expect(electronCount({ atoms: parseXyz(preset('oh').xyz), charge: -1, multiplicity: 1 })).toBe(10);
  });
});

describe('presets', () => {
  it('have consistent charge and multiplicity', () => {
    for (const p of PRESETS) {
      const n = electronCount({ atoms: parseXyz(p.xyz), charge: p.charge, multiplicity: p.multiplicity });
      expect((n + p.multiplicity - 1) % 2, p.id).toBe(0);
      expect(p.bases.length, p.id).toBeGreaterThan(0);
    }
    expect(() => preset('nope')).toThrow();
  });

  it('converge with Hartree-Fock in their largest basis', () => {
    const session = new qchem.Session();
    for (const p of PRESETS) {
      const basis = p.bases[p.bases.length - 1];
      const out = session.compute({
        atoms: parseXyz(p.xyz),
        charge: p.charge,
        multiplicity: p.multiplicity,
        basis,
        method: 'hf',
        warmStart: false,
      }) as ComputeOutput;
      expect(out.converged, p.id).toBe(true);
      expect(out.energy, p.id).toBeLessThan(0);
    }
    session.free();
  });

  it('give the known STO-3G energy of benzene', () => {
    const out = compute({ atoms: parseXyz(preset('c6h6').xyz), charge: 0, multiplicity: 1, basis: 'sto-3g', method: 'hf' });
    expect(out.nbf).toBe(36);
    // PySCF RHF/STO-3G at r(CC) = 1.397 Å, r(CH) = 1.084 Å.
    expect(out.energy).toBeCloseTo(-227.89, 1);
  });
});

/** Integrate f over a cube [-L, L]^3 with the midpoint rule. */
function integrate(f: (x: number, y: number, z: number) => number, L: number, n: number): number {
  const h = (2 * L) / n;
  let s = 0;
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) s += f(-L + (i + 0.5) * h, -L + (j + 0.5) * h, -L + (k + 0.5) * h);
  return s * h * h * h;
}

describe('orbital evaluation (CPU reference for the GPU shader)', () => {
  const water = compute(inputFor(refCase('h2o', '6-31g*')));
  const oh = compute(inputFor(refCase('oh', '6-31g')));

  it('integrates the density to the electron count', () => {
    const h2 = compute({ atoms: parseXyz(preset('h2').xyz), charge: 0, multiplicity: 1, basis: '6-31g', method: 'hf' });
    expect(integrate((x, y, z) => densityAt(h2, x, y, z), 7, 70)).toBeCloseTo(2, 3);
    // A cube grid cannot resolve the tight oxygen 1s functions, so heavy atoms only get a loose check.
    const ao = new Float64Array(water.nbf);
    expect(Math.abs(integrate((x, y, z) => densityAt(water, x, y, z, ao), 7, 90) - 10)).toBeLessThan(0.15);
    expect(Math.abs(integrate((x, y, z) => densityAt(oh, x, y, z), 7, 90) - 9)).toBeLessThan(0.15);
  });

  it('gives normalized molecular orbitals', () => {
    const ao = new Float64Array(water.nbf);
    const lumo = water.nalpha;
    const norm = integrate(
      (x, y, z) => {
        evalAo(water.basis, x, y, z, ao);
        return moValue(ao, water.ca, water.nmo, lumo) ** 2;
      },
      8,
      70,
    );
    expect(norm).toBeCloseTo(1, 2);
  });

  it('uses the PySCF cartesian component order', () => {
    // A d shell's xy component must vanish on the x axis and the xx component must not.
    const ao = new Float64Array(water.nbf);
    const shells = water.basis.shells;
    const d = Array.from({ length: shells.length / 8 }, (_, s) => s).find((s) => shells[8 * s + 3] === 2)!;
    const [cx, cy, cz] = shells.slice(8 * d, 8 * d + 3);
    const f = shells[8 * d + 6];
    evalAo(water.basis, cx + 0.5, cy, cz, ao);
    expect(ao[f]).toBeGreaterThan(0);
    expect(ao[f + 1]).toBe(0);
    expect(ao[f + 3]).toBe(0);
    evalAo(water.basis, cx, cy + 0.4, cz + 0.4, ao);
    // Each cartesian component is unit-normalized, so yz is sqrt(3) larger than yy here.
    expect(ao[f + 4]).toBeCloseTo(ao[f + 3] * Math.sqrt(3), 12);
  });

  it('labels orbital levels and frontier orbitals', () => {
    const levels = orbitalLevels(water);
    expect(levels[4]).toMatchObject({ label: 'HOMO', occupied: true });
    expect(levels[5]).toMatchObject({ label: 'LUMO', occupied: false });
    expect(levels[3].label).toBe('HOMO−1');
    expect(levels[6].label).toBe('LUMO+1');
    const { homo, lumo } = frontier(water);
    expect(homo).toBe(water.epsA[4]);
    expect(lumo).toBe(water.epsA[5]);
    const open = frontier(oh);
    expect(open.homo).toBe(Math.max(oh.epsA[4], oh.epsB[3]));
    expect(orbitalLevels(oh, 'beta')[3].label).toBe('HOMO');
  });
});
