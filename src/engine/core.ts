import { PARALLEL_MIN_ATOMS, type Densities, type EriPart } from './eri-pool';
import type { ComputeInput, ComputeOutput, WorkerRequest, WorkerResponse } from './types';

/** The parts of the WASM module the engine uses, so tests can pass the real module loaded with `initSync`. */
export interface QchemModule {
  Session: new () => {
    compute(input: unknown, eriParts?: unknown[]): unknown;
    densities(): unknown;
    gradient(twoElectron: Float64Array): unknown;
    reset(): void;
  };
  Optimizer: new (natom: number) => {
    step(coords: Float64Array, energy: number, gradient: Float64Array): Float64Array | undefined;
    free(): void;
  };
}

const yieldToEvents = () => new Promise((r) => setTimeout(r, 0));

/** Computes the two-electron integrals in parallel (see `EriPool`). */
export interface EriSource {
  readonly size: number;
  parts(input: ComputeInput): Promise<EriPart[]>;
  gradient(input: ComputeInput, densities: Densities): Promise<Float64Array>;
}

/** Request handler shared by the worker and the tests. */
export function createEngine(qchem: QchemModule, post: (msg: WorkerResponse) => void, eris: EriSource | null = null) {
  const session = new qchem.Session();
  let activeRelax = 0;

  const parallel = (input: ComputeInput) => eris !== null && eris.size > 1 && input.atoms.length >= PARALLEL_MIN_ATOMS;

  // One calculation at a time: with the integrals computed elsewhere, a second request could
  // otherwise run on the session while the first waits, and the gradient would mix them up.
  // Small molecules run synchronously, as they always have.
  let lock: Promise<unknown> = Promise.resolve();
  const compute = (input: ComputeInput): ComputeOutput | Promise<ComputeOutput> => {
    if (!parallel(input)) return session.compute(input) as ComputeOutput;
    const run = lock.then(() => computeParallel(input, eris!));
    lock = run.catch(() => undefined);
    return run;
  };

  const computeParallel = async (input: ComputeInput, eris: EriSource): Promise<ComputeOutput> => {
    const t0 = performance.now();
    const out = session.compute({ ...input, gradient: false }, await eris.parts(input)) as ComputeOutput;
    // The time includes the integrals computed by the other workers.
    if (!input.gradient) return { ...out, millis: performance.now() - t0 };
    const twoElectron = await eris.gradient(input, session.densities() as Densities);
    const g = session.gradient(twoElectron) as { gradient: number[]; maxForce: number };
    return { ...out, gradient: g.gradient, maxForce: g.maxForce, millis: performance.now() - t0 };
  };

  async function relax(id: number, input: ComputeInput, maxSteps: number) {
    activeRelax = id;
    const opt = new qchem.Optimizer(input.atoms.length);
    let atoms = input.atoms;
    try {
      for (let step = 0; step < maxSteps; step++) {
        const pending = compute({ ...input, atoms, gradient: true, history: false });
        const out = pending instanceof Promise ? await pending : pending;
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
          post({ id: req.id, type: 'result', output: await compute(req.input) });
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
        case 'configure':
          // Handled by the worker, which owns the integral workers.
          break;
      }
    } catch (e) {
      post({ id: req.id, type: 'error', message: e instanceof Error ? e.message : String(e) });
    }
  };
}
