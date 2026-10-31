import { electronCount, parseXyz, toXyz, type Atom } from './chem/molecule';
import { PRESETS, preset } from './chem/presets';
import { orbitalLevels } from './chem/orbitals';
import { EngineClient } from './engine/client';
import type { BasisName, ComputeInput, ComputeOutput, MethodName } from './engine/types';
import { OrbitCamera } from './render/camera';
import { add, dot, raySphere, scale, sub, type Vec3 } from './render/math';
import { Renderer, type ViewSettings } from './render/renderer';
import { volumeInputs, type OrbitalSelection } from './render/volume-data';
import { readouts } from './ui/format';
import { describeOrbital, legendElements } from './ui/legend';
import { renderDiagram } from './ui/orbital-diagram';
import { atomRadius, boundingSphere, buildScene } from './ui/scene';

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

const ui = {
  canvas: $<HTMLCanvasElement>('view'),
  noGpu: $('no-gpu'),
  preset: $<HTMLSelectElement>('preset'),
  xyzBox: $<HTMLDetailsElement>('xyz-box'),
  xyz: $<HTMLTextAreaElement>('xyz'),
  loadXyz: $<HTMLButtonElement>('load-xyz'),
  charge: $<HTMLInputElement>('charge'),
  multiplicity: $<HTMLInputElement>('multiplicity'),
  method: $<HTMLSelectElement>('method'),
  basis: $<HTMLSelectElement>('basis'),
  relax: $<HTMLButtonElement>('relax'),
  resetGeometry: $<HTMLButtonElement>('reset-geometry'),
  showCloud: $<HTMLInputElement>('show-cloud'),
  cloudOpacity: $<HTMLInputElement>('cloud-opacity'),
  cloudMin: $<HTMLInputElement>('cloud-min'),
  showOrbital: $<HTMLInputElement>('show-orbital'),
  iso: $<HTMLInputElement>('iso'),
  status: $('status'),
  readouts: $('readouts'),
  orbitalName: $('orbital-name'),
  diagram: document.getElementById('diagram') as unknown as SVGSVGElement,
  saveImage: $<HTMLButtonElement>('save-image'),
  playScf: $<HTMLButtonElement>('play-scf'),
  scfFrame: $<HTMLInputElement>('scf-frame'),
  scfLabel: $('scf-label'),
  help: $<HTMLButtonElement>('help'),
  intro: $<HTMLDialogElement>('intro'),
  legendCloud: $('legend-cloud'),
  legendOrbital: $('legend-orbital'),
  legendOrbitalName: $('legend-orbital-name'),
  legendAtoms: $('legend-atoms'),
  legendHelp: $<HTMLButtonElement>('legend-help'),
};

interface State {
  atoms: Atom[];
  /** Geometry to return to with Reset. */
  initialAtoms: Atom[];
  charge: number;
  multiplicity: number;
  basis: BasisName;
  method: MethodName;
  result: ComputeOutput | null;
  orbital: OrbitalSelection | null;
  /** SCF iteration shown, or null for the converged result. */
  frame: number | null;
  relaxId: number | null;
  hover: number | null;
}

const state: State = {
  atoms: [],
  initialAtoms: [],
  charge: 0,
  multiplicity: 1,
  basis: '6-31g',
  method: 'hf',
  result: null,
  orbital: null,
  frame: null,
  relaxId: null,
  hover: null,
};

const camera = new OrbitCamera();
const engine = new EngineClient(new Worker(new URL('./engine/engine.worker.ts', import.meta.url), { type: 'module' }));
let renderer: Renderer | null = null;
/** `?headless` draws offscreen only; headless browsers cannot present a WebGPU canvas. */
const headless = new URLSearchParams(location.search).has('headless');
/** Grid points along the longest box edge; `?res=64` lowers it for slow GPUs. */
const volumeResolution = Math.max(16, Math.min(160, Number(new URLSearchParams(location.search).get('res')) || 112));
let needsRender = true;
let playTimer: number | null = null;

function viewSettings(): ViewSettings {
  return {
    showCloud: ui.showCloud.checked,
    showOrbital: ui.showOrbital.checked && state.orbital !== null,
    cloudOpacity: 40 * Number(ui.cloudOpacity.value) ** 2,
    logMin: Number(ui.cloudMin.value),
    logMax: 1,
    iso: Number(ui.iso.value),
    orbitalOpacity: 0.55,
  };
}

function setStatus(text: string, kind: 'ok' | 'busy' | 'error' = 'ok') {
  ui.status.textContent = text;
  ui.status.className = `status ${kind === 'ok' ? '' : kind}`;
}

function input(extra: Partial<ComputeInput> = {}): ComputeInput {
  return {
    atoms: state.atoms.map((a) => ({ symbol: a.symbol, xyz: a.xyz })),
    charge: state.charge,
    multiplicity: state.multiplicity,
    basis: state.basis,
    method: state.method,
    ...extra,
  };
}

// ---- scene and volume -------------------------------------------------------------------------

function updateScene() {
  if (!renderer) return;
  const { spheres, bonds } = buildScene(state.atoms, state.hover);
  renderer.setAtoms(spheres, bonds);
  needsRender = true;
}

function updateVolume() {
  const res = state.result;
  if (!renderer || !res) return;
  const it = state.frame === null ? null : res.history[state.frame];
  try {
    renderer.updateVolume(volumeInputs(res, state.orbital, it ?? res, volumeResolution));
  } catch (e) {
    setStatus(e instanceof Error ? e.message : String(e), 'error');
  }
  needsRender = true;
}

function frameLoop() {
  if (renderer && needsRender) {
    needsRender = false;
    resizeCanvas();
    renderer.render(camera, viewSettings());
  }
  requestAnimationFrame(frameLoop);
}

function resizeCanvas() {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(1, Math.floor(ui.canvas.clientWidth * dpr));
  const h = Math.max(1, Math.floor(ui.canvas.clientHeight * dpr));
  if (ui.canvas.width !== w || ui.canvas.height !== h) {
    ui.canvas.width = w;
    ui.canvas.height = h;
  }
}

// ---- results ----------------------------------------------------------------------------------

function defaultOrbital(res: ComputeOutput): OrbitalSelection {
  return { spin: 'alpha', index: res.nalpha - 1 };
}

function showResult(res: ComputeOutput, keepHistory = false) {
  const sameShape = state.result && state.result.nmo === res.nmo && state.result.unrestricted === res.unrestricted;
  state.result = res;
  if (!sameShape || !state.orbital || state.orbital.index >= res.nmo) state.orbital = defaultOrbital(res);
  if (!keepHistory) stopPlayback();
  state.frame = null;
  ui.scfFrame.max = String(Math.max(0, res.history.length - 1));
  ui.scfFrame.value = ui.scfFrame.max;
  ui.scfFrame.disabled = res.history.length === 0;
  ui.scfLabel.textContent = res.history.length ? `iteration ${res.history.length} of ${res.history.length}` : 'SCF iterations';
  ui.readouts.replaceChildren(
    ...readouts(res).flatMap((r) => {
      const dt = document.createElement('dt');
      dt.textContent = r.label;
      const dd = document.createElement('dd');
      dd.textContent = r.value;
      if (r.title) dd.title = r.title;
      return [dt, dd];
    }),
  );
  drawDiagram();
  updateVolume();
}

function drawDiagram() {
  const res = state.result;
  if (!res) return;
  renderDiagram(ui.diagram, res, state.orbital, (sel) => {
    state.orbital = sel;
    ui.showOrbital.checked = true;
    drawDiagram();
    updateVolume();
  });
  if (state.orbital) {
    const level = orbitalLevels(res, state.orbital.spin)[state.orbital.index];
    const name = `${level.label}${res.unrestricted ? (state.orbital.spin === 'alpha' ? ' α' : ' β') : ''}`;
    ui.orbitalName.textContent = `· ${name}`;
    ui.legendOrbitalName.textContent = `${name} orbital`;
    ui.legendOrbital.title = `${name}: ${describeOrbital(level.label)}. Blue where its wavefunction is positive, orange where it is negative; where they meet is a node, where this orbital's electrons are never found.`;
  }
  updateLegend();
}

function updateLegend() {
  ui.legendCloud.hidden = !ui.showCloud.checked;
  ui.legendOrbital.hidden = !(ui.showOrbital.checked && state.orbital);
  ui.legendAtoms.replaceChildren(
    ...legendElements(state.atoms).map((e) => {
      const item = document.createElement('span');
      item.className = 'legend-item';
      item.title = e.name;
      const swatch = document.createElement('span');
      swatch.className = 'swatch';
      swatch.style.background = e.color;
      item.append(swatch, e.symbol);
      return item;
    }),
  );
}

function describe(res: ComputeOutput) {
  return `${res.energy.toFixed(6)} Ha · ${res.iterations} iterations · ${Math.round(res.millis)} ms`;
}

async function recompute(opts: { history?: boolean } = {}) {
  if (state.atoms.length === 0) return;
  setStatus('Solving the SCF…', 'busy');
  try {
    const res = await engine.computeLatest(input({ history: opts.history ?? true }));
    if (!res) return;
    showResult(res);
    setStatus(res.converged ? describe(res) : `SCF did not converge (${res.iterations} iterations)`, res.converged ? 'ok' : 'error');
  } catch (e) {
    setStatus(e instanceof Error ? e.message : String(e), 'error');
  }
}

// ---- molecule setup ---------------------------------------------------------------------------

function setMolecule(atoms: Atom[], charge: number, multiplicity: number) {
  cancelRelax();
  state.atoms = atoms;
  state.initialAtoms = atoms.map((a) => ({ symbol: a.symbol, xyz: [...a.xyz] as Vec3 }));
  state.charge = charge;
  state.multiplicity = multiplicity;
  state.result = null;
  state.orbital = null;
  ui.charge.value = String(charge);
  ui.multiplicity.value = String(multiplicity);
  ui.xyz.value = toXyz(atoms).split('\n').slice(2).join('\n');
  syncMethodOptions();
  engine.reset();
  const { center, radius } = boundingSphere(atoms);
  camera.fit(center, radius);
  updateScene();
  updateLegend();
  recompute();
}

function loadPreset(id: string) {
  const p = preset(id);
  for (const opt of Array.from(ui.basis.options)) opt.disabled = !p.bases.includes(opt.value as BasisName);
  if (!p.bases.includes(state.basis)) {
    state.basis = p.bases[p.bases.length - 1];
    ui.basis.value = state.basis;
  }
  setMolecule(parseXyz(p.xyz), p.charge, p.multiplicity);
}

/** Open shells only work with UHF for now. */
function syncMethodOptions() {
  const open = state.multiplicity !== 1;
  for (const opt of Array.from(ui.method.options)) opt.disabled = open && opt.value !== 'hf';
  if (open && state.method !== 'hf') {
    state.method = 'hf';
    ui.method.value = 'hf';
  }
}

// ---- relax ------------------------------------------------------------------------------------

function cancelRelax() {
  if (state.relaxId !== null) engine.cancelRelax(state.relaxId);
  state.relaxId = null;
  ui.relax.textContent = 'Relax geometry';
}

function relax() {
  if (state.relaxId !== null) return cancelRelax();
  ui.relax.textContent = 'Stop';
  setStatus('Relaxing…', 'busy');
  state.relaxId = engine.relax(input({ history: false }), {
    onStep: (res, step) => {
      state.atoms = state.atoms.map((a, i) => ({ symbol: a.symbol, xyz: res.coordinates.slice(3 * i, 3 * i + 3) as Vec3 }));
      updateScene();
      showResult(res);
      setStatus(`Relaxing: step ${step + 1}, E = ${res.energy.toFixed(6)} Ha, max force ${res.maxForce?.toExponential(1)}`, 'busy');
    },
    onDone: (converged, steps) => {
      state.relaxId = null;
      ui.relax.textContent = 'Relax geometry';
      setStatus(converged ? `Relaxed in ${steps} steps · ${describe(state.result!)}` : `Stopped after ${steps} steps without converging`, converged ? 'ok' : 'error');
    },
    onError: (message) => {
      state.relaxId = null;
      ui.relax.textContent = 'Relax geometry';
      setStatus(message, 'error');
    },
  });
}

// ---- SCF playback -----------------------------------------------------------------------------

function stopPlayback() {
  if (playTimer !== null) clearInterval(playTimer);
  playTimer = null;
  ui.playScf.textContent = '▶ Play SCF';
}

function showFrame(k: number) {
  const res = state.result;
  if (!res || res.history.length === 0) return;
  state.frame = k >= res.history.length - 1 ? null : k;
  ui.scfFrame.value = String(k);
  const h = res.history[k];
  ui.scfLabel.textContent = `iteration ${k + 1} of ${res.history.length} · E = ${h.energy.toFixed(6)} · error ${h.error.toExponential(1)}`;
  updateVolume();
}

async function playScf() {
  if (playTimer !== null) return stopPlayback();
  cancelRelax();
  // Replay from the atomic-density guess, which is more instructive than a warm start.
  setStatus('Recording SCF iterations…', 'busy');
  const res = await engine.compute(input({ history: true, warmStart: false })).catch((e: Error) => {
    setStatus(e.message, 'error');
    return null;
  });
  if (!res) return;
  showResult(res, true);
  setStatus(describe(res));
  let k = 0;
  ui.playScf.textContent = '■ Stop';
  showFrame(0);
  playTimer = window.setInterval(() => {
    k++;
    if (k >= res.history.length) return stopPlayback();
    showFrame(k);
  }, 450);
}

// ---- pointer interaction ----------------------------------------------------------------------

function ndc(ev: PointerEvent | WheelEvent): [number, number] {
  const r = ui.canvas.getBoundingClientRect();
  return [((ev.clientX - r.left) / r.width) * 2 - 1, 1 - ((ev.clientY - r.top) / r.height) * 2];
}

function pick(ev: PointerEvent): number | null {
  const { origin, dir } = camera.ray(...ndc(ev));
  let best: number | null = null;
  let bestT = Infinity;
  state.atoms.forEach((a, i) => {
    const t = raySphere(origin, dir, a.xyz, atomRadius(a.symbol) * 1.1);
    if (t !== null && t < bestT) {
      bestT = t;
      best = i;
    }
  });
  return best;
}

let drag: { kind: 'orbit'; x: number; y: number } | { kind: 'atom'; index: number; planePoint: Vec3; normal: Vec3 } | null = null;

ui.canvas.addEventListener('pointerdown', (ev) => {
  ui.canvas.setPointerCapture(ev.pointerId);
  const hit = pick(ev);
  if (hit !== null) {
    cancelRelax();
    stopPlayback();
    drag = { kind: 'atom', index: hit, planePoint: [...state.atoms[hit].xyz] as Vec3, normal: camera.basis().forward };
  } else drag = { kind: 'orbit', x: ev.clientX, y: ev.clientY };
  ui.canvas.classList.add('dragging');
});

ui.canvas.addEventListener('pointermove', (ev) => {
  if (!drag) {
    const hover = pick(ev);
    ui.canvas.classList.toggle('over-atom', hover !== null);
    if (hover !== state.hover) {
      state.hover = hover;
      updateScene();
    }
    return;
  }
  if (drag.kind === 'orbit') {
    camera.orbit(ev.clientX - drag.x, ev.clientY - drag.y);
    drag.x = ev.clientX;
    drag.y = ev.clientY;
    needsRender = true;
    return;
  }
  const { origin, dir } = camera.ray(...ndc(ev));
  const denom = dot(dir, drag.normal);
  if (Math.abs(denom) < 1e-6) return;
  const t = dot(sub(drag.planePoint, origin), drag.normal) / denom;
  const p = add(origin, scale(dir, t));
  const index = drag.index;
  state.atoms = state.atoms.map((a, i) => (i === index ? { symbol: a.symbol, xyz: p } : a));
  updateScene();
  recompute({ history: false });
});

const endDrag = () => {
  drag = null;
  ui.canvas.classList.remove('dragging');
};
ui.canvas.addEventListener('pointerup', endDrag);
ui.canvas.addEventListener('pointercancel', endDrag);
ui.canvas.addEventListener(
  'wheel',
  (ev) => {
    ev.preventDefault();
    camera.zoom(Math.exp(ev.deltaY * 0.001));
    needsRender = true;
  },
  { passive: false },
);
window.addEventListener('resize', () => (needsRender = true));

// ---- controls ---------------------------------------------------------------------------------

for (const p of PRESETS) ui.preset.add(new Option(p.name, p.id));
ui.preset.add(new Option('Custom (pasted XYZ)', 'custom'));
ui.preset.addEventListener('change', () => {
  if (ui.preset.value === 'custom') ui.xyzBox.open = true;
  else loadPreset(ui.preset.value);
});
ui.loadXyz.addEventListener('click', () => {
  try {
    const atoms = parseXyz(ui.xyz.value);
    if (atoms.length === 0) throw new Error('No atoms found. Use lines like "O 0 0 0" in Angstrom.');
    ui.preset.value = 'custom';
    for (const opt of Array.from(ui.basis.options)) opt.disabled = false;
    const n = electronCount({ atoms, charge: state.charge, multiplicity: 1 });
    setMolecule(atoms, state.charge, n % 2 === 0 ? 1 : 2);
  } catch (e) {
    setStatus(e instanceof Error ? e.message : String(e), 'error');
  }
});
ui.charge.addEventListener('change', () => {
  state.charge = Math.round(Number(ui.charge.value));
  const n = electronCount({ atoms: state.atoms, charge: state.charge, multiplicity: 1 });
  // Keep the spin state valid for the new electron count.
  if ((n + state.multiplicity - 1) % 2 !== 0) state.multiplicity = state.multiplicity === 1 ? 2 : state.multiplicity - 1;
  ui.multiplicity.value = String(state.multiplicity);
  syncMethodOptions();
  recompute();
});
ui.multiplicity.addEventListener('change', () => {
  state.multiplicity = Math.max(1, Math.round(Number(ui.multiplicity.value)));
  syncMethodOptions();
  recompute();
});
ui.method.addEventListener('change', () => {
  state.method = ui.method.value as MethodName;
  recompute();
});
ui.basis.addEventListener('change', () => {
  state.basis = ui.basis.value as BasisName;
  recompute();
});
ui.relax.addEventListener('click', relax);
ui.resetGeometry.addEventListener('click', () => {
  cancelRelax();
  state.atoms = state.initialAtoms.map((a) => ({ symbol: a.symbol, xyz: [...a.xyz] as Vec3 }));
  updateScene();
  recompute();
});
for (const el of [ui.cloudOpacity, ui.cloudMin, ui.iso]) el.addEventListener('input', () => (needsRender = true));
ui.showCloud.addEventListener('input', () => {
  updateLegend();
  needsRender = true;
});
ui.showOrbital.addEventListener('input', () => {
  if (ui.showOrbital.checked && !state.orbital && state.result) state.orbital = defaultOrbital(state.result);
  updateLegend();
  needsRender = true;
});
ui.playScf.addEventListener('click', playScf);
ui.saveImage.addEventListener('click', async () => {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const a = document.createElement('a');
  a.href = await snapshot(Math.round(ui.canvas.clientWidth * dpr), Math.round(ui.canvas.clientHeight * dpr));
  a.download = `${ui.preset.value}-${state.method}-${state.basis}.png`;
  a.click();
});
ui.scfFrame.addEventListener('input', () => {
  stopPlayback();
  showFrame(Number(ui.scfFrame.value));
});

// ---- intro ------------------------------------------------------------------------------------

const INTRO_SEEN = 'quantum-molecules:intro-seen';

function showIntro(section?: string) {
  if (!ui.intro.open) ui.intro.showModal();
  if (section) document.getElementById(section)?.scrollIntoView({ block: 'start' });
  else ui.intro.scrollTop = 0;
}

ui.help.addEventListener('click', () => showIntro());
ui.legendHelp.addEventListener('click', () => showIntro('intro-colors'));
ui.intro.addEventListener('close', () => {
  try {
    localStorage.setItem(INTRO_SEEN, '1');
  } catch {
    // Storage can be blocked; the intro then just shows again next time.
  }
});
// A click on the backdrop lands on the dialog element itself.
ui.intro.addEventListener('click', (ev) => {
  if (ev.target === ui.intro) ui.intro.close();
});

function introSeen(): boolean {
  try {
    return localStorage.getItem(INTRO_SEEN) !== null;
  } catch {
    return false;
  }
}

// ---- start ------------------------------------------------------------------------------------

async function start() {
  // Tests run with ?headless, which skips the intro unless ?intro is also given.
  if ((!headless || new URLSearchParams(location.search).has('intro')) && !introSeen()) showIntro();
  ui.method.value = state.method;
  ui.basis.value = state.basis;
  try {
    renderer = await Renderer.create(headless ? null : ui.canvas);
    renderer.device.lost.then((info) => setStatus(`The GPU device was lost: ${info.message}`, 'error'));
  } catch (e) {
    ui.noGpu.hidden = false;
    ui.noGpu.textContent = `${e instanceof Error ? e.message : String(e)} The calculations still run and the numbers on the right are real, but the 3D view needs a browser with WebGPU (recent Chrome, Edge, or Safari).`;
  }
  ui.preset.value = 'h2o';
  await engine.ready;
  loadPreset('h2o');
  requestAnimationFrame(frameLoop);
}

/** Render the current view to a PNG data URL. */
async function snapshot(width = ui.canvas.clientWidth || 800, height = ui.canvas.clientHeight || 600): Promise<string> {
  if (!renderer) throw new Error('WebGPU is not available.');
  camera.aspect = width / height;
  const pixels = await renderer.snapshot(camera, viewSettings(), width, height);
  const c = document.createElement('canvas');
  c.width = width;
  c.height = height;
  c.getContext('2d')!.putImageData(new ImageData(pixels, width, height), 0, 0);
  return c.toDataURL('image/png');
}

// Exposed for end-to-end tests and the console.
Object.assign(window, {
  __qm: {
    state,
    engine,
    camera,
    snapshot,
    volumeResolution,
    get renderer() {
      return renderer;
    },
  },
});

start();
