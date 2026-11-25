/// <reference lib="webworker" />
import init, * as qchem from '@qchem/qchem_wasm.js';
import { createEngine, type EriSource } from './core';
import { EriPool } from './eri-pool';
import type { WorkerRequest, WorkerResponse } from './types';

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const post = (msg: WorkerResponse) => ctx.postMessage(msg);

// The integral workers, started when the page says how many it wants (`configure`). This worker
// only waits while they run.
let pool: EriPool | null = null;
const eris: EriSource = {
  get size() {
    return pool?.size ?? 1;
  },
  parts: (input) => pool!.parts(input),
  gradient: (input, densities) => pool!.gradient(input, densities),
};

const engine = init().then(() => {
  post({ id: 0, type: 'ready' });
  return createEngine(qchem, post, eris);
});

ctx.onmessage = async (ev: MessageEvent<WorkerRequest>) => {
  const req = ev.data;
  if (req.type === 'configure') {
    if (!pool && req.threads > 1)
      pool = new EriPool(req.threads, () => new Worker(new URL('./eri.worker.ts', import.meta.url), { type: 'module' }));
    return;
  }
  const handle = await engine;
  handle(req);
};
