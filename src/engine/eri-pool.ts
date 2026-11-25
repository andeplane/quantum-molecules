import type { ComputeInput } from './types';

/** One stripe of the two-electron integrals, as returned by the WASM `eriPart`. */
export interface EriPart {
  index: Uint32Array;
  value: Float64Array;
}

/** Alpha and beta density matrices (nbf x nbf), from `Session.densities`. */
export interface Densities {
  da: Float64Array;
  db: Float64Array;
}

export type EriPartRequest =
  | { id: number; kind: 'integrals'; input: ComputeInput; part: number; stride: number }
  | { id: number; kind: 'gradient'; input: ComputeInput; densities: Densities; part: number; stride: number };

export type EriPartResponse =
  | ({ id: number; kind: 'integrals' } & EriPart)
  | { id: number; kind: 'gradient'; gradient: Float64Array }
  | { id: number; error: string };

/**
 * Below this many atoms the integrals take a few milliseconds and handing them out costs more
 * than it saves.
 */
export const PARALLEL_MIN_ATOMS = 8;

/**
 * Spreads the two-electron integrals, the slowest step of every calculation, over several
 * workers. Each worker computes every `stride`-th row of shell pairs, and the engine assembles
 * the stripes into the full tensor. Nothing else is shared, so the result is identical to a
 * single-threaded run. The two-electron part of the nuclear gradient is shared out the same way.
 */
export class EriPool {
  private workers: Worker[];
  private nextId = 1;
  private pending = new Map<number, { resolve: (msg: EriPartResponse) => void; reject: (e: Error) => void }>();

  constructor(size: number, spawn: () => Worker) {
    this.workers = Array.from({ length: size }, () => {
      const w = spawn();
      w.onmessage = (ev: MessageEvent<EriPartResponse>) => {
        const msg = ev.data;
        const p = this.pending.get(msg.id);
        if (!p) return;
        this.pending.delete(msg.id);
        if ('error' in msg) p.reject(new Error(msg.error));
        else p.resolve(msg);
      };
      return w;
    });
  }

  get size() {
    return this.workers.length;
  }

  private run(make: (part: number, stride: number) => Omit<EriPartRequest, 'id'>): Promise<EriPartResponse[]> {
    const stride = this.workers.length;
    return Promise.all(
      this.workers.map(
        (w, part) =>
          new Promise<EriPartResponse>((resolve, reject) => {
            const id = this.nextId++;
            this.pending.set(id, { resolve, reject });
            w.postMessage({ id, ...make(part, stride) } as EriPartRequest);
          }),
      ),
    );
  }

  /** Every stripe of the two-electron integrals. */
  async parts(input: ComputeInput): Promise<EriPart[]> {
    const out = await this.run((part, stride) => ({ kind: 'integrals', input, part, stride }));
    return out.map((m) => {
      if (!('index' in m)) throw new Error('unexpected reply from an integral worker');
      return { index: m.index, value: m.value };
    });
  }

  /** The two-electron part of the nuclear gradient (bohr), summed over the workers. */
  async gradient(input: ComputeInput, densities: Densities): Promise<Float64Array> {
    const out = await this.run((part, stride) => ({ kind: 'gradient', input, densities, part, stride }));
    const sum = new Float64Array(3 * input.atoms.length);
    for (const m of out) {
      if (!('gradient' in m)) throw new Error('unexpected reply from an integral worker');
      m.gradient.forEach((v, i) => (sum[i] += v));
    }
    return sum;
  }
}
