import { describe, expect, it } from 'vitest';
import { createEngine } from '../src/engine/core';
import { parseXyz } from '../src/chem/molecule';
import type { Densities, EriPart } from '../src/engine/eri-pool';
import type { ComputeInput, WorkerResponse } from '../src/engine/types';
import { inputFor, qchem, refCase } from './helpers/qchem';

function engine() {
  const messages: WorkerResponse[] = [];
  const handle = createEngine(qchem, (m) => messages.push(m));
  return { handle, messages };
}

describe('engine request handler', () => {
  it('answers compute requests with results', async () => {
    const { handle, messages } = engine();
    const c = refCase('h2o', 'sto-3g');
    await handle({ id: 7, type: 'compute', input: inputFor(c) });
    expect(messages).toHaveLength(1);
    const msg = messages[0];
    expect(msg.type).toBe('result');
    if (msg.type === 'result') expect(msg.output.energy).toBeCloseTo(c.energy, 8);
    expect(msg.id).toBe(7);
  });

  it('turns engine errors into error messages', async () => {
    const { handle, messages } = engine();
    await handle({ id: 3, type: 'compute', input: { ...inputFor(refCase('h2o', 'sto-3g')), multiplicity: 2 } });
    expect(messages).toEqual([{ id: 3, type: 'error', message: expect.stringMatching(/multiplicity/) }]);
  });

  it('streams relaxation steps until converged', async () => {
    const { handle, messages } = engine();
    const input = { ...inputFor(refCase('h2o_distorted', '6-31g*')), basis: 'sto-3g' as const };
    await handle({ id: 1, type: 'relax', input });
    const steps = messages.filter((m) => m.type === 'relax-step');
    const done = messages[messages.length - 1];
    expect(done).toMatchObject({ type: 'relax-done', converged: true, steps: steps.length });
    const energies = steps.map((m) => (m.type === 'relax-step' ? m.output.energy : 0));
    expect(energies[energies.length - 1]).toBeLessThan(energies[0]);
    expect(Math.min(...energies)).toBeCloseTo(energies[energies.length - 1], 6);
    for (const m of steps) if (m.type === 'relax-step') expect(m.output.gradient).toHaveLength(9);
  });

  it('stops after maxSteps without converging', async () => {
    const { handle, messages } = engine();
    const input = { ...inputFor(refCase('h2o_distorted', '6-31g*')), basis: 'sto-3g' as const };
    await handle({ id: 1, type: 'relax', input, maxSteps: 2 });
    expect(messages.map((m) => m.type)).toEqual(['relax-step', 'relax-step', 'relax-done']);
    expect(messages[2]).toMatchObject({ converged: false, steps: 2 });
  });

  it('cancels a running relaxation between steps', async () => {
    const { handle, messages } = engine();
    const input = { ...inputFor(refCase('h2o_distorted', '6-31g*')), basis: 'sto-3g' as const };
    const running = handle({ id: 5, type: 'relax', input });
    await handle({ id: 6, type: 'cancel' });
    await running;
    expect(messages.filter((m) => m.type === 'relax-step')).toHaveLength(1);
    expect(messages.some((m) => m.type === 'relax-done')).toBe(false);
  });

  it('resets the warm-start orbitals', async () => {
    const { handle, messages } = engine();
    const input = inputFor(refCase('h2o', '6-31g*'));
    await handle({ id: 1, type: 'compute', input });
    await handle({ id: 2, type: 'compute', input });
    await handle({ id: 3, type: 'reset' });
    await handle({ id: 4, type: 'compute', input });
    const iters = messages.map((m) => (m.type === 'result' ? m.output.iterations : 0));
    expect(iters[1]).toBeLessThan(iters[0]);
    expect(iters[2]).toBe(iters[0]);
  });
});

describe('parallel integrals', () => {
  // Stands in for the worker pool: the same WASM calls, run in turn.
  const stripes = 3;
  const pool = {
    size: stripes,
    parts: async (input: ComputeInput) => Array.from({ length: stripes }, (_, k) => qchem.eriPart(input, k, stripes) as EriPart),
    gradient: async (input: ComputeInput, d: Densities) => {
      const sum = new Float64Array(3 * input.atoms.length);
      for (let k = 0; k < stripes; k++) qchem.eriGradientPart(input, d.da, d.db, k, stripes).forEach((v, i) => (sum[i] += v));
      return sum;
    },
  };
  // Eleven atoms, enough to take the parallel path.
  const propane: ComputeInput = {
    atoms: parseXyz(
      'C 0 0 0\nC 1.27 0.85 0\nC -1.27 0.85 0\nH 0 -0.64 0.88\nH 0 -0.64 -0.88\nH 1.27 1.50 0.88\nH 1.27 1.50 -0.88\nH 2.16 0.21 0\nH -1.27 1.50 0.88\nH -1.27 1.50 -0.88\nH -2.16 0.21 0',
    ),
    charge: 0,
    multiplicity: 1,
    basis: 'sto-3g',
    method: 'hf',
    gradient: true,
  };

  it('gives the same energy and forces as one worker', async () => {
    const serial = engine();
    const messages: WorkerResponse[] = [];
    const parallel = createEngine(qchem, (m) => messages.push(m), pool);
    await serial.handle({ id: 1, type: 'compute', input: propane });
    await parallel({ id: 1, type: 'compute', input: propane });
    const a = serial.messages[0];
    const b = messages[0];
    if (a.type !== 'result' || b.type !== 'result') throw new Error('expected results');
    expect(b.output.energy).toBeCloseTo(a.output.energy, 10);
    b.output.gradient!.forEach((g, i) => expect(g).toBeCloseTo(a.output.gradient![i], 9));
  });
});
