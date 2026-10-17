/// <reference lib="webworker" />
import init, * as qchem from '@qchem/qchem_wasm.js';
import { createEngine } from './core';
import type { WorkerRequest, WorkerResponse } from './types';

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const post = (msg: WorkerResponse) => ctx.postMessage(msg);
const engine = init().then(() => {
  post({ id: 0, type: 'ready' });
  return createEngine(qchem, post);
});

ctx.onmessage = async (ev: MessageEvent<WorkerRequest>) => {
  const handle = await engine;
  handle(ev.data);
};
