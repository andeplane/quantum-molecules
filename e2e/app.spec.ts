import { expect, test, type Page } from '@playwright/test';
import type { ComputeOutput } from '../src/engine/types';
import { evalPacked, volumeInputs } from '../src/render/volume-data';

declare global {
  interface Window {
    __qm: {
      state: { result: ComputeOutput | null; atoms: { symbol: string; xyz: [number, number, number] }[]; orbital: { spin: string; index: number } | null; frame: number | null; relaxId: number | null };
      snapshot(width?: number, height?: number): Promise<string>;
      volumeResolution: number;
      renderer: { readVolume(): Promise<{ dims: [number, number, number]; data: Float32Array } | null> } | null;
      camera: { viewProjection: Float32Array };
    };
  }
}

async function open(page: Page, query = '') {
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(e.message));
  page.on('console', (m) => {
    if (m.type() === 'error' && !m.text().includes('404')) errors.push(m.text());
  });
  await page.goto(`/?headless${query}`);
  await page.waitForFunction(() => window.__qm?.state.result && window.__qm.renderer, null, { timeout: 60_000 });
  return errors;
}

/** Drag an atom sideways on screen by `dx` pixels. */
async function dragAtom(page: Page, index: number, dx = 40) {
  const pos = await page.evaluate((index) => {
    const m = window.__qm.camera.viewProjection;
    const p = window.__qm.state.atoms[index].xyz;
    const x = m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12];
    const y = m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13];
    const w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
    const r = document.getElementById('view')!.getBoundingClientRect();
    return { x: r.left + ((x / w + 1) / 2) * r.width, y: r.top + ((1 - y / w) / 2) * r.height };
  }, index);
  await page.mouse.move(pos.x, pos.y);
  await page.mouse.down();
  for (let k = 1; k <= 5; k++) await page.mouse.move(pos.x + (dx / 5) * k, pos.y);
  await page.mouse.up();
}

/** The orbital-diagram line labelled HOMO or LUMO. */
const level = (page: Page, label: string) => page.locator('#diagram .level').filter({ has: page.locator('text.label', { hasText: new RegExp(`^${label}$`) }) });

const openSettings = (page: Page) => page.click('#settings summary');

const result = (page: Page) => page.evaluate(() => window.__qm.state.result);
const current = async (page: Page) => (await result(page))!;
const energyText = (page: Page) => page.locator('#readouts dd').first().textContent();

/** Fraction of pixels that differ from the background colour. */
async function coverage(page: Page): Promise<number> {
  return page.evaluate(async () => {
    const url = await window.__qm.snapshot(320, 240);
    const img = new Image();
    img.src = url;
    await img.decode();
    const c = document.createElement('canvas');
    c.width = 320;
    c.height = 240;
    const ctx = c.getContext('2d')!;
    ctx.drawImage(img, 0, 0);
    const d = ctx.getImageData(0, 0, 320, 240).data;
    let n = 0;
    for (let i = 0; i < d.length; i += 4) if (d[i] + d[i + 1] + d[i + 2] > 60) n++;
    return n / (320 * 240);
  });
}

test('computes water on load and draws it', async ({ page }) => {
  const errors = await open(page);
  const res = await current(page);
  // RHF/6-31G water at the preset geometry.
  expect(res.converged).toBe(true);
  expect(res.energy).toBeCloseTo(-75.98397, 4);
  await expect(page.locator('#energy')).toHaveText('-75.98397 Ha');
  await expect(page.locator('#status')).toContainText('Solved in');
  // Element labels sit over the atoms.
  await expect(page.locator('#labels span')).toHaveText(['O', 'H', 'H']);
  expect(await energyText(page)).toContain('-75.98397');
  expect(await page.locator('#diagram .level').count()).toBeGreaterThan(5);
  expect(await coverage(page)).toBeGreaterThan(0.05);
  expect(errors).toEqual([]);
});

test('the GPU volume matches the CPU reference', async ({ page }) => {
  await open(page, '&res=40');
  const res = await current(page);
  const orbital = await page.evaluate(() => window.__qm.state.orbital);
  const vol = await page.evaluate(async () => {
    const v = await window.__qm.renderer!.readVolume();
    return v && { dims: v.dims, data: Array.from(v.data) };
  });
  expect(vol).not.toBeNull();
  const inputs = volumeInputs(res, orbital as never, res, 40);
  expect(vol!.dims).toEqual(inputs.box.dims);
  const [nx, ny, nz] = inputs.box.dims;
  let checked = 0;
  for (let z = 0; z < nz; z += 3)
    for (let y = 0; y < ny; y += 3)
      for (let x = 0; x < nx; x += 3) {
        const p = [0, 1, 2].map((k) => inputs.box.min[k] + ([x, y, z][k] * (inputs.box.max[k] - inputs.box.min[k])) / (inputs.box.dims[k] - 1)) as [number, number, number];
        const cpu = evalPacked(inputs, p);
        const i = 2 * ((z * ny + y) * nx + x);
        // Half-float storage: about 3 significant digits.
        expect(Math.abs(vol!.data[i] - cpu.rho)).toBeLessThan(2e-3 * Math.abs(cpu.rho) + 1e-6);
        expect(Math.abs(vol!.data[i + 1] - cpu.psi)).toBeLessThan(2e-3 * Math.abs(cpu.psi) + 1e-6);
        checked++;
      }
  expect(checked).toBeGreaterThan(1000);
});

test('switching molecule, method and basis recomputes', async ({ page }) => {
  await open(page);
  await openSettings(page);
  await page.selectOption('#preset', 'ch4');
  await expect.poll(async () => (await result(page))?.nbf).toBe(17);
  await page.selectOption('#method', 'pbe');
  await expect.poll(async () => (await result(page))?.energy).toBeLessThan(-40.4);
  await page.selectOption('#basis', 'sto-3g');
  await expect.poll(async () => (await result(page))?.nbf).toBe(9);
  // Open shells only allow Hartree-Fock: choosing O2 switches the method back.
  await page.selectOption('#preset', 'o2');
  await expect(page.locator('#method')).toHaveValue('hf');
  await expect.poll(async () => (await result(page))?.unrestricted).toBe(true);
  expect(await page.locator('#method option[value=pbe]').evaluate((o) => (o as HTMLOptionElement).disabled)).toBe(true);
  // Benzene only offers STO-3G.
  await page.selectOption('#preset', 'c6h6');
  await expect(page.locator('#basis')).toHaveValue('sto-3g');
  await expect.poll(async () => (await result(page))?.nbf, { timeout: 30_000 }).toBe(36);
});

test('relaxes the geometry', async ({ page }) => {
  await open(page);
  await openSettings(page);
  await page.selectOption('#basis', 'sto-3g');
  await expect.poll(async () => (await result(page))?.nbf).toBe(7);
  const before = (await current(page)).energy;
  await page.click('#relax');
  await expect(page.locator('#status')).toContainText('Relaxed in', { timeout: 60_000 });
  const after = await current(page);
  expect(after.energy).toBeLessThan(before);
  expect(after.maxForce!).toBeLessThan(4.5e-4);
  await expect(page.locator('#relax')).toHaveText('Relax geometry');
});

test('dragging an atom moves it and re-runs the SCF', async ({ page }) => {
  await open(page);
  const before = await current(page);
  await dragAtom(page, 1);
  await expect.poll(async () => (await result(page))?.coordinates[3 * 1]).not.toBeCloseTo(before.coordinates[3], 3);
  const after = await current(page);
  expect(after.energy).not.toBeCloseTo(before.energy, 5);
  expect(after.energy).toBeGreaterThan(before.energy);
  await page.click('#reset-geometry');
  await expect.poll(async () => (await result(page))?.energy).toBeCloseTo(before.energy, 6);
});

test('selecting an orbital level and playing the SCF', async ({ page }) => {
  await open(page);
  // The view starts on the cloud alone; clicking a level shows its orbital.
  await expect(page.locator('#orbital-name')).toHaveText('');
  await expect(page.locator('#show-orbital')).not.toBeChecked();
  await level(page, 'HOMO').click();
  await expect(page.locator('#orbital-name')).toContainText('HOMO');
  await expect(page.locator('#show-orbital')).toBeChecked();
  await page.locator('#diagram .level:not(.occupied)').first().click();
  await expect(page.locator('#orbital-name')).toContainText('LUMO');
  await page.click('#play-scf');
  await expect(page.locator('#scf-label')).toContainText('step 1 of');
  await expect(page.locator('#play-scf')).toHaveText('▶ Watch it solve', { timeout: 20_000 });
  await page.locator('#scf-frame').fill('0');
  await expect.poll(() => page.evaluate(() => window.__qm.state.frame)).toBe(0);
  expect(await coverage(page)).toBeGreaterThan(0.05);
});

test('loads pasted XYZ and reports errors', async ({ page }) => {
  await open(page);
  await openSettings(page);
  await page.click('#xyz-box summary');
  await page.fill('#xyz', 'H 0 0 0\nF 0 0 0.92');
  await page.click('#load-xyz');
  await expect.poll(async () => (await result(page))?.nbf).toBe(11);
  await expect(page.locator('#preset')).toHaveValue('custom');
  await page.fill('#xyz', 'Na 0 0 0');
  await page.click('#load-xyz');
  await expect(page.locator('#status')).toHaveClass(/error/);
});

test('explains itself on the first visit and on demand', async ({ page }) => {
  await open(page, '&intro');
  const intro = page.locator('#intro');
  await expect(intro).toBeVisible();
  await expect(intro).toContainText('where the electrons are');
  await expect(page.locator('#tip')).toBeHidden();
  await page.locator('#intro-close').click();
  await expect(intro).toBeHidden();
  await expect.poll(() => page.evaluate(() => localStorage.getItem('quantum-molecules:intro-seen'))).toBe('1');

  // It can be reopened, and closed with Escape.
  await page.locator('#help').click();
  await expect(intro).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(intro).toBeHidden();
  // "What does this mean?" in the legend opens the detailed part.
  await page.locator('#legend-help').click();
  await expect(page.locator('#intro-more')).toHaveAttribute('open', '');
  await page.keyboard.press('Escape');

  // The legend names what is on screen and follows the view.
  await expect(page.locator('#legend-title')).toHaveText('Water H₂O');
  await expect(page.locator('#legend-atoms')).toHaveText('OH');
  await expect(page.locator('#legend-orbital')).toBeHidden();
  await level(page, 'LUMO').click();
  await expect(page.locator('#legend-orbital-name')).toHaveText('LUMO orbital');
  await page.locator('#show-orbital').uncheck();
  await expect(page.locator('#legend-orbital')).toBeHidden();
  await page.selectOption('#preset', 'ch4');
  await expect(page.locator('#legend-atoms')).toHaveText('CH');
  await expect(page.locator('#preset-blurb')).toContainText('tetrahedron');
});

test('tips walk a newcomer through the page', async ({ page }) => {
  await open(page, '&intro');
  await page.locator('#intro-close').click();
  const tip = page.locator('#tip-text');
  await expect(tip).toContainText('drag one of the atoms');
  await dragAtom(page, 1);
  await expect(tip).toContainText('Relax geometry');
  await page.click('#relax');
  await expect(tip).toContainText('Energy levels', { timeout: 60_000 });
  await level(page, 'HOMO').click();
  await expect(tip).toContainText('another molecule');
  await page.selectOption('#preset', 'nh3');
  await expect(page.locator('#tip')).toBeHidden();
  await expect.poll(() => page.evaluate(() => localStorage.getItem('quantum-molecules:tips-done'))).toBe('1');

  // They can be restarted from the help dialog, and hidden.
  await page.locator('#help').click();
  await page.locator('#restart-tips').click();
  await expect(page.locator('#intro')).toBeHidden();
  await expect(tip).toContainText('drag one of the atoms');
  await page.locator('#tip-close').click();
  await expect(page.locator('#tip')).toBeHidden();
});

test('does not show the intro or tips to returning visitors', async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem('quantum-molecules:intro-seen', '1');
    localStorage.setItem('quantum-molecules:tips-done', '1');
  });
  await open(page, '&intro');
  await expect(page.locator('#intro')).toBeHidden();
  await expect(page.locator('#tip')).toBeHidden();
});
