import { describe, expect, it } from 'vitest';
import { REFERENCE, compute, inputFor, qchem, refCase } from './helpers/qchem';

const maxDiff = (a: ArrayLike<number>, b: ArrayLike<number>) => {
  expect(a.length).toBe(b.length);
  let m = 0;
  for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
};

describe('WASM engine against PySCF', () => {
  // The DFT references use a much finer grid than the interactive default, and the interactive
  // SCF convergence (1e-5 in the orbital gradient) limits first-order properties to about 1e-5.
  for (const c of REFERENCE) {
    const dft = c.xc !== null;
    it(`${c.name} ${c.basis} ${c.xc ?? c.method}`, () => {
      const out = compute(inputFor(c, { gradient: true, grid: 'ultrafine' }));
      expect(out.converged).toBe(true);
      expect(out.nbf).toBe(c.nbf);
      expect(out.unrestricted).toBe(c.method === 'uhf');
      expect(Math.abs(out.energy - c.energy)).toBeLessThan(dft ? 2e-5 : 1e-8);
      expect(maxDiff(out.dipole, c.dipole)).toBeLessThan(1e-4);
      expect(maxDiff(out.gradient!, c.gradient.flat())).toBeLessThan(1e-4);
      expect(out.maxForce).toBeCloseTo(Math.max(...c.gradient.flat().map(Math.abs)), 4);
      if (c.moEnergies) expect(maxDiff(out.epsA, c.moEnergies)).toBeLessThan(1e-4);
    });
  }
});

describe('Session', () => {
  const water = refCase('h2o', '6-31g*');

  it('returns orbitals in the documented layout', () => {
    const out = compute(inputFor(water));
    expect(out.ca.length).toBe(out.nbf * out.nmo);
    expect(out.epsA.length).toBe(out.nmo);
    expect(out.nalpha).toBe(5);
    expect(out.nbeta).toBe(5);
    expect(out.coordinates).toEqual(water.atoms.flatMap((a) => a.xyz));
    expect(out.gradient).toBeUndefined();
    expect(out.maxForce).toBeUndefined();
    expect(out.history).toEqual([]);
    expect(out.millis).toBeGreaterThanOrEqual(0);
  });

  it('warm-starts from the previous geometry', () => {
    const session = new qchem.Session();
    const cold = compute(inputFor(water), session);
    const moved = water.atoms.map((a, i) => ({ ...a, xyz: [a.xyz[0], a.xyz[1] + (i === 1 ? 0.02 : 0), a.xyz[2]] as [number, number, number] }));
    const warm = compute({ ...inputFor(water), atoms: moved }, session);
    const fresh = compute({ ...inputFor(water), atoms: moved, warmStart: false }, session);
    expect(warm.converged).toBe(true);
    expect(warm.iterations).toBeLessThan(fresh.iterations);
    expect(Math.abs(warm.energy - fresh.energy)).toBeLessThan(1e-8);
    expect(cold.energy).toBeLessThan(warm.energy);
    session.free();
  });

  it('ignores stale orbitals after switching basis, molecule or spin', () => {
    const session = new qchem.Session();
    compute(inputFor(water), session);
    const sto = compute(inputFor(refCase('h2o', 'sto-3g')), session);
    expect(Math.abs(sto.energy - refCase('h2o', 'sto-3g').energy)).toBeLessThan(1e-8);
    const n2 = refCase('n2', '6-31g*');
    expect(Math.abs(compute(inputFor(n2), session).energy - n2.energy)).toBeLessThan(1e-8);
    const oh = refCase('oh', '6-31g');
    expect(Math.abs(compute(inputFor(oh), session).energy - oh.energy)).toBeLessThan(1e-8);
    session.reset();
    session.free();
  });

  it('keeps orbital signs stable while an atom is dragged', () => {
    // The HOMO is drawn with blue and orange lobes; a sign flip between frames swaps them.
    const session = new qchem.Session();
    const input = { ...inputFor(refCase('h2o', 'sto-3g')), basis: '6-31g' as const };
    let prev = compute(input, session);
    let atoms = input.atoms.map((a) => ({ ...a, xyz: [...a.xyz] as [number, number, number] }));
    for (let step = 0; step < 10; step++) {
      atoms = atoms.map((a, i) => (i === 1 ? { ...a, xyz: [a.xyz[0] + 0.01, a.xyz[1] + 0.04, a.xyz[2]] } : a));
      const next = compute({ ...input, atoms }, session);
      const homo = next.nalpha - 1;
      // Same basis-function layout, small move: coefficient vectors should point the same way.
      for (const k of [homo, homo + 1]) {
        let d = 0;
        for (let i = 0; i < next.nbf; i++) d += prev.ca[i * next.nmo + k] * next.ca[i * next.nmo + k];
        expect(d, `step ${step}, orbital ${k}`).toBeGreaterThan(0);
      }
      prev = next;
    }
    session.free();
  });

  it('records the SCF history for playback', () => {
    const out = compute(inputFor(water, { history: true }));
    expect(out.history.length).toBe(out.iterations);
    const last = out.history[out.history.length - 1];
    expect(Math.abs(last.energy - out.energy)).toBeLessThan(1e-6);
    expect(last.ca.length).toBe(out.ca.length);
    expect(out.history[0].error).toBeGreaterThan(last.error);
  });

  it('exposes the basis in the GPU layout', () => {
    const out = compute(inputFor(water));
    const { shells, prims, componentNorms } = out.basis;
    expect(shells.length % 8).toBe(0);
    const nshell = shells.length / 8;
    // 6-31G* water: O has 1s, 2sp, 3sp, d; each H has 2 s shells.
    expect(nshell).toBe(3 + 2 + 1 + 4);
    const last = 8 * (nshell - 1);
    expect(shells[last + 6] + [1, 3, 6][shells[last + 3]]).toBe(out.nbf);
    expect(prims.length / 2).toBe(shells.reduce((n, _, i) => (i % 8 === 4 ? n + shells[i] : n), 0));
    expect(componentNorms.map((n) => n.length)).toEqual([1, 3, 6]);
  });

  it('reports bad input as errors', () => {
    const session = new qchem.Session();
    const base = inputFor(water);
    expect(() => session.compute({ ...base, atoms: [{ symbol: 'Xx', xyz: [0, 0, 0] }] })).toThrow(/Xx/);
    expect(() => session.compute({ ...base, basis: 'cc-pvdz' })).toThrow(/basis/i);
    expect(() => session.compute({ ...base, method: 'b3lyp' })).toThrow(/method|functional/i);
    expect(() => session.compute({ ...base, multiplicity: 2 })).toThrow();
    expect(() => session.compute({ ...base, multiplicity: 3, method: 'pbe' })).toThrow(/open-shell|restricted|unrestricted/i);
    expect(() => session.compute({ atoms: [] })).toThrow();
    session.free();
  });

  it('knows which elements each basis covers', () => {
    expect(qchem.hasBasis('sto-3g', 'Ne')).toBe(true);
    expect(qchem.hasBasis('6-31g*', 'C')).toBe(true);
    expect(qchem.hasBasis('6-31g', 'Na')).toBe(false);
    expect(qchem.hasBasis('nope', 'H')).toBe(false);
  });
});

describe('Optimizer', () => {
  it('relaxes distorted water to the STO-3G minimum', () => {
    const session = new qchem.Session();
    let atoms = refCase('h2o_distorted', '6-31g*').atoms;
    const opt = new qchem.Optimizer(atoms.length);
    let converged = false;
    let energy = 0;
    for (let step = 0; step < 50 && !converged; step++) {
      const out = compute({ atoms, charge: 0, multiplicity: 1, basis: 'sto-3g', method: 'hf', gradient: true }, session);
      energy = out.energy;
      const next = opt.step(Float64Array.from(out.coordinates), out.energy, Float64Array.from(out.gradient!));
      if (!next) converged = true;
      else atoms = atoms.map((a, i) => ({ symbol: a.symbol, xyz: [next[3 * i], next[3 * i + 1], next[3 * i + 2]] }));
    }
    expect(converged).toBe(true);
    expect(opt.steps).toBeLessThan(30);
    const [o, h1, h2] = atoms.map((a) => a.xyz);
    const sub = (a: number[], b: number[]) => a.map((x, k) => x - b[k]);
    const len = (v: number[]) => Math.hypot(...v);
    const r1 = len(sub(h1, o)) * 0.52917721092;
    const r2 = len(sub(h2, o)) * 0.52917721092;
    const cos = sub(h1, o).reduce((s, x, k) => s + x * sub(h2, o)[k], 0) / (len(sub(h1, o)) * len(sub(h2, o)));
    expect(r1).toBeCloseTo(0.989, 2);
    expect(r2).toBeCloseTo(0.989, 2);
    expect((Math.acos(cos) * 180) / Math.PI).toBeCloseTo(100.0, 0);
    expect(energy).toBeCloseTo(-74.9659, 3);
    opt.free();
    session.free();
  });
});
