import type { ComputeInput, ComputeOutput, WorkerRequest, WorkerResponse } from './types';

export interface RelaxCallbacks {
  onStep(output: ComputeOutput, step: number): void;
  onDone(converged: boolean, steps: number): void;
  onError(message: string): void;
}

/** Talks to the engine worker. Interactive computes are "latest wins": while one runs, only the newest pending request is kept. */
export class EngineClient {
  private worker: Worker;
  private nextId = 1;
  private pending = new Map<number, { resolve: (o: ComputeOutput) => void; reject: (e: Error) => void }>();
  private relaxes = new Map<number, RelaxCallbacks>();
  private busy = false;
  private queued: { input: ComputeInput; resolve: (o: ComputeOutput | null) => void } | null = null;
  readonly ready: Promise<void>;

  constructor(worker: Worker) {
    this.worker = worker;
    let markReady: () => void;
    this.ready = new Promise((r) => (markReady = r));
    worker.onmessage = (ev: MessageEvent<WorkerResponse>) => {
      const msg = ev.data;
      if (msg.type === 'ready') return markReady();
      if (msg.type === 'relax-step' || msg.type === 'relax-done') {
        const cb = this.relaxes.get(msg.id);
        if (!cb) return;
        if (msg.type === 'relax-step') cb.onStep(msg.output, msg.step);
        else {
          this.relaxes.delete(msg.id);
          cb.onDone(msg.converged, msg.steps);
        }
        return;
      }
      if (msg.type === 'error' && this.relaxes.has(msg.id)) {
        const cb = this.relaxes.get(msg.id)!;
        this.relaxes.delete(msg.id);
        cb.onError(msg.message);
        return;
      }
      const p = this.pending.get(msg.id);
      if (!p) return;
      this.pending.delete(msg.id);
      if (msg.type === 'result') p.resolve(msg.output);
      else if (msg.type === 'error') p.reject(new Error(msg.message));
    };
  }

  private send(req: WorkerRequest) {
    this.worker.postMessage(req);
  }

  compute(input: ComputeInput): Promise<ComputeOutput> {
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.send({ id, type: 'compute', input });
    });
  }

  /**
   * Compute, but if a computation is already running keep only the latest request.
   * Superseded requests resolve with null.
   */
  computeLatest(input: ComputeInput): Promise<ComputeOutput | null> {
    if (this.busy) {
      this.queued?.resolve(null);
      return new Promise((resolve) => (this.queued = { input, resolve }));
    }
    this.busy = true;
    return this.compute(input).finally(() => {
      this.busy = false;
      const q = this.queued;
      this.queued = null;
      if (q) this.computeLatest(q.input).then(q.resolve, () => q.resolve(null));
    });
  }

  relax(input: ComputeInput, callbacks: RelaxCallbacks, maxSteps = 60): number {
    const id = this.nextId++;
    this.relaxes.set(id, callbacks);
    this.send({ id, type: 'relax', input, maxSteps });
    return id;
  }

  cancelRelax(id: number) {
    this.relaxes.delete(id);
    this.send({ id, type: 'cancel' });
  }

  reset() {
    this.send({ id: this.nextId++, type: 'reset' });
  }
}
