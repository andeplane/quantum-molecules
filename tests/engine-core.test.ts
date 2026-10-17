import { describe, expect, it } from 'vitest';
import { createEngine } from '../src/engine/core';
import type { WorkerResponse } from '../src/engine/types';
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
