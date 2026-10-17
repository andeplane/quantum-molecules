import type { ComputeInput, ComputeOutput, WorkerRequest, WorkerResponse } from './types';

/** The parts of the WASM module the engine uses, so tests can pass the real module loaded with `initSync`. */
export interface QchemModule {
  Session: new () => { compute(input: unknown): unknown; reset(): void };
  Optimizer: new (natom: number) => {
    step(coords: Float64Array, energy: number, gradient: Float64Array): Float64Array | undefined;
    free(): void;
  };
}

const yieldToEvents = () => new Promise((r) => setTimeout(r, 0));

/** Request handler shared by the worker and the tests. */
export function createEngine(qchem: QchemModule, post: (msg: WorkerResponse) => void) {
  const session = new qchem.Session();
  let activeRelax = 0;

  const compute = (input: ComputeInput) => session.compute(input) as ComputeOutput;

  async function relax(id: number, input: ComputeInput, maxSteps: number) {
    activeRelax = id;
    const opt = new qchem.Optimizer(input.atoms.length);
    let atoms = input.atoms;
    try {
      for (let step = 0; step < maxSteps; step++) {
        const out = compute({ ...input, atoms, gradient: true, history: false });
        if (activeRelax !== id) return;
        post({ id, type: 'relax-step', output: out, step });
        const next = opt.step(Float64Array.from(out.coordinates), out.energy, Float64Array.from(out.gradient ?? []));
        if (!next) {
          post({ id, type: 'relax-done', converged: true, steps: step + 1 });
          return;
        }
        atoms = atoms.map((a, i) => ({ symbol: a.symbol, xyz: [next[3 * i], next[3 * i + 1], next[3 * i + 2]] }));
        // Let a cancel message in between steps.
        await yieldToEvents();
        if (activeRelax !== id) return;
      }
      post({ id, type: 'relax-done', converged: false, steps: maxSteps });
    } finally {
      opt.free();
      if (activeRelax === id) activeRelax = 0;
    }
  }

  return async function handle(req: WorkerRequest): Promise<void> {
    try {
      switch (req.type) {
        case 'compute':
          post({ id: req.id, type: 'result', output: compute(req.input) });
          break;
        case 'relax':
          await relax(req.id, req.input, req.maxSteps ?? 60);
          break;
        case 'cancel':
          activeRelax = 0;
          break;
        case 'reset':
          session.reset();
          break;
      }
    } catch (e) {
      post({ id: req.id, type: 'error', message: e instanceof Error ? e.message : String(e) });
    }
  };
}
