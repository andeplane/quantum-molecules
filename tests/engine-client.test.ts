import { describe, expect, it } from 'vitest';
import { EngineClient } from '../src/engine/client';
import type { ComputeInput, ComputeOutput, WorkerRequest, WorkerResponse } from '../src/engine/types';

/** A stand-in for the worker that records requests and lets the test reply. */
class FakeWorker {
  onmessage: ((ev: MessageEvent<WorkerResponse>) => void) | null = null;
  sent: WorkerRequest[] = [];
  postMessage(req: WorkerRequest) {
    this.sent.push(req);
  }
  reply(msg: WorkerResponse) {
    this.onmessage?.({ data: msg } as MessageEvent<WorkerResponse>);
  }
}

const input = (charge: number): ComputeInput => ({ atoms: [], charge, multiplicity: 1, basis: 'sto-3g', method: 'hf' });
const output = (energy: number) => ({ energy }) as ComputeOutput;
const flush = () => new Promise((r) => setTimeout(r, 0));

function client() {
  const worker = new FakeWorker();
  return { worker, client: new EngineClient(worker as unknown as Worker) };
}

describe('EngineClient', () => {
  it('resolves ready when the worker reports it', async () => {
    const { worker, client: c } = client();
    worker.reply({ id: 0, type: 'ready' });
    await expect(c.ready).resolves.toBeUndefined();
  });

  it('matches results and errors to requests by id', async () => {
    const { worker, client: c } = client();
    const a = c.compute(input(0));
    const b = c.compute(input(1));
    const [ra, rb] = worker.sent;
    worker.reply({ id: rb.id, type: 'error', message: 'nope' });
    worker.reply({ id: ra.id, type: 'result', output: output(-1) });
    await expect(a).resolves.toEqual(output(-1));
    await expect(b).rejects.toThrow('nope');
  });

  it('keeps only the newest request while a computation runs', async () => {
    const { worker, client: c } = client();
    const first = c.computeLatest(input(0));
    const second = c.computeLatest(input(1));
    const third = c.computeLatest(input(2));
    expect(worker.sent).toHaveLength(1);
    await expect(second).resolves.toBeNull();
    worker.reply({ id: worker.sent[0].id, type: 'result', output: output(0) });
    await expect(first).resolves.toEqual(output(0));
    await flush();
    expect(worker.sent).toHaveLength(2);
    const req = worker.sent[1];
    expect(req.type === 'compute' && req.input.charge).toBe(2);
    worker.reply({ id: req.id, type: 'result', output: output(2) });
    await expect(third).resolves.toEqual(output(2));
  });

  it('keeps going after a failed latest-wins computation', async () => {
    const { worker, client: c } = client();
    const first = c.computeLatest(input(0));
    worker.reply({ id: worker.sent[0].id, type: 'error', message: 'bad' });
    await expect(first).rejects.toThrow('bad');
    const next = c.computeLatest(input(1));
    expect(worker.sent).toHaveLength(2);
    worker.reply({ id: worker.sent[1].id, type: 'result', output: output(1) });
    await expect(next).resolves.toEqual(output(1));
  });

  it('routes relaxation steps and completion to callbacks', () => {
    const { worker, client: c } = client();
    const steps: number[] = [];
    const done: boolean[] = [];
    const id = c.relax(input(0), { onStep: (_, s) => steps.push(s), onDone: (ok) => done.push(ok), onError: () => {} }, 10);
    expect(worker.sent[0]).toMatchObject({ id, type: 'relax', maxSteps: 10 });
    worker.reply({ id, type: 'relax-step', output: output(0), step: 0 });
    worker.reply({ id, type: 'relax-step', output: output(0), step: 1 });
    worker.reply({ id, type: 'relax-done', converged: true, steps: 2 });
    worker.reply({ id, type: 'relax-step', output: output(0), step: 2 });
    expect(steps).toEqual([0, 1]);
    expect(done).toEqual([true]);
  });

  it('reports relaxation errors and ignores steps after cancel', () => {
    const { worker, client: c } = client();
    const errors: string[] = [];
    const steps: number[] = [];
    const cb = { onStep: (_: ComputeOutput, s: number) => steps.push(s), onDone: () => {}, onError: (m: string) => errors.push(m) };
    const a = c.relax(input(0), cb);
    worker.reply({ id: a, type: 'error', message: 'scf failed' });
    expect(errors).toEqual(['scf failed']);
    const b = c.relax(input(0), cb);
    c.cancelRelax(b);
    expect(worker.sent[worker.sent.length - 1]).toEqual({ id: b, type: 'cancel' });
    worker.reply({ id: b, type: 'relax-step', output: output(0), step: 0 });
    expect(steps).toEqual([]);
  });

  it('sends reset requests', () => {
    const { worker, client: c } = client();
    c.reset();
    expect(worker.sent[0].type).toBe('reset');
  });
});
