/// <reference lib="webworker" />
import init, { eriGradientPart, eriPart } from '@qchem/qchem_wasm.js';
import type { EriPart, EriPartRequest, EriPartResponse } from './eri-pool';

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const ready = init();

ctx.onmessage = async (ev: MessageEvent<EriPartRequest>) => {
  await ready;
  const req = ev.data;
  try {
    if (req.kind === 'integrals') {
      const out = eriPart(req.input, req.part, req.stride) as EriPart;
      const msg: EriPartResponse = { id: req.id, kind: 'integrals', index: out.index, value: out.value };
      ctx.postMessage(msg, [out.index.buffer, out.value.buffer]);
    } else {
      const { da, db } = req.densities;
      const gradient = eriGradientPart(req.input, da, db, req.part, req.stride);
      const msg: EriPartResponse = { id: req.id, kind: 'gradient', gradient };
      ctx.postMessage(msg, [gradient.buffer]);
    }
  } catch (e) {
    const msg: EriPartResponse = { id: req.id, error: e instanceof Error ? e.message : String(e) };
    ctx.postMessage(msg);
  }
};
