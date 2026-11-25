/// <reference lib="webworker" />
import init, * as qchem from '@qchem/qchem_wasm.js';
import { createEngine } from './core';
import { EriPool } from './eri-pool';
import type { WorkerRequest, WorkerResponse } from './types';

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const post = (msg: WorkerResponse) => ctx.postMessage(msg);
// One integral worker per core; this worker only waits while they run. `?threads=1` turns it off.
const threads = Number(new URLSearchParams(location.search).get('threads')) || navigator.hardwareConcurrency || 4;
const poolSize = Math.max(1, Math.min(16, threads));
const pool =
  poolSize > 1 ? new EriPool(poolSize, () => new Worker(new URL('./eri.worker.ts', import.meta.url), { type: 'module' })) : null;
const engine = init().then(() => {
  post({ id: 0, type: 'ready' });
  return createEngine(qchem, post, pool);
});

ctx.onmessage = async (ev: MessageEvent<WorkerRequest>) => {
  const handle = await engine;
  handle(ev.data);
};
