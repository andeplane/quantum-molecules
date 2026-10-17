import type { ComputeOutput } from '../engine/types';

export interface OrbitalSelection {
  spin: 'alpha' | 'beta';
  index: number;
}

export interface VolumeBox {
  min: [number, number, number];
  max: [number, number, number];
  dims: [number, number, number];
}

/** Buffers for `volume.wgsl`, packed exactly as the shader reads them. */
export interface VolumeInputs {
  nbf: number;
  nshell: number;
  nocc: number;
  /** 32 bytes per shell: center (3 f32), l, nprim, prim, fn0 (u32), minExp (f32). */
  shells: ArrayBuffer;
  /** exponent, coefficient per primitive. */
  prims: Float32Array<ArrayBuffer>;
  norms: Float32Array<ArrayBuffer>;
  /** (nocc + 1) rows of nbf coefficients; the last row is the selected orbital. */
  coefs: Float32Array<ArrayBuffer>;
  weights: Float32Array<ArrayBuffer>;
  box: VolumeBox;
}

export const MAX_BF = 128;
/** Padding around the nuclei, in bohr. Diffuse virtual orbitals reach further than the density. */
export const BOX_PADDING = 5;

export function volumeBox(coordinates: ArrayLike<number>, maxDim = 128, padding = BOX_PADDING): VolumeBox {
  const min: [number, number, number] = [Infinity, Infinity, Infinity];
  const max: [number, number, number] = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < coordinates.length; i += 3)
    for (let k = 0; k < 3; k++) {
      min[k] = Math.min(min[k], coordinates[i + k] - padding);
      max[k] = Math.max(max[k], coordinates[i + k] + padding);
    }
  const extent = max.map((m, k) => m - min[k]);
  const h = Math.max(...extent) / (maxDim - 1);
  const dims = extent.map((e) => Math.min(maxDim, Math.ceil(e / h - 1e-9) + 1)) as [number, number, number];
  // Grow the box so the spacing is the same along every axis.
  for (let k = 0; k < 3; k++) {
    const c = (min[k] + max[k]) / 2;
    const half = ((dims[k] - 1) * h) / 2;
    min[k] = c - half;
    max[k] = c + half;
  }
  return { min, max, dims };
}

export function volumeInputs(
  res: ComputeOutput,
  orbital: OrbitalSelection | null,
  coefficients: { ca: ArrayLike<number>; cb: ArrayLike<number> } = res,
  maxDim = 128,
): VolumeInputs {
  const { nbf, nmo } = res;
  if (nbf > MAX_BF) throw new Error(`The renderer supports at most ${MAX_BF} basis functions (this system has ${nbf}).`);
  const { shells: sh, prims: pr, componentNorms } = res.basis;
  const nshell = sh.length / 8;
  const shells = new ArrayBuffer(32 * nshell);
  const f32 = new Float32Array(shells);
  const u32 = new Uint32Array(shells);
  for (let s = 0; s < nshell; s++) {
    const o = 8 * s;
    const nprim = sh[o + 4];
    const p0 = sh[o + 5];
    let minExp = Infinity;
    for (let p = 0; p < nprim; p++) minExp = Math.min(minExp, pr[2 * (p0 + p)]);
    f32[o] = sh[o];
    f32[o + 1] = sh[o + 1];
    f32[o + 2] = sh[o + 2];
    u32[o + 3] = sh[o + 3];
    u32[o + 4] = nprim;
    u32[o + 5] = p0;
    u32[o + 6] = sh[o + 6];
    f32[o + 7] = minExp;
  }
  const cb = res.unrestricted ? coefficients.cb : coefficients.ca;
  const rows: { c: ArrayLike<number>; k: number; w: number }[] = [];
  for (let k = 0; k < res.nalpha; k++) rows.push({ c: coefficients.ca, k, w: res.unrestricted ? 1 : 2 });
  if (res.unrestricted) for (let k = 0; k < res.nbeta; k++) rows.push({ c: cb, k, w: 1 });
  const nocc = rows.length;
  const coefs = new Float32Array((nocc + 1) * nbf);
  rows.forEach(({ c, k }, r) => {
    for (let i = 0; i < nbf; i++) coefs[r * nbf + i] = c[i * nmo + k];
  });
  if (orbital) {
    const c = orbital.spin === 'beta' ? cb : coefficients.ca;
    for (let i = 0; i < nbf; i++) coefs[nocc * nbf + i] = c[i * nmo + orbital.index];
  }
  return {
    nbf,
    nshell,
    nocc,
    shells,
    prims: Float32Array.from(pr),
    norms: Float32Array.from(componentNorms.flat()),
    coefs,
    weights: Float32Array.from(rows.map((r) => r.w)),
    box: volumeBox(res.coordinates, maxDim),
  };
}

/** What one invocation of `volume.wgsl` computes, run on the CPU from the packed buffers. */
export function evalPacked(v: VolumeInputs, p: [number, number, number]): { rho: number; psi: number } {
  const f32 = new Float32Array(v.shells);
  const u32 = new Uint32Array(v.shells);
  const ao = new Float64Array(v.nbf);
  const n = v.norms;
  for (let s = 0; s < v.nshell; s++) {
    const o = 8 * s;
    const dx = p[0] - f32[o];
    const dy = p[1] - f32[o + 1];
    const dz = p[2] - f32[o + 2];
    const r2 = dx * dx + dy * dy + dz * dz;
    const l = u32[o + 3];
    const f = u32[o + 6];
    let radial = 0;
    if (f32[o + 7] * r2 < 40) for (let k = 0; k < u32[o + 4]; k++) radial += v.prims[2 * (u32[o + 5] + k) + 1] * Math.exp(-v.prims[2 * (u32[o + 5] + k)] * r2);
    if (l === 0) ao[f] = radial * n[0];
    else if (l === 1) [dx, dy, dz].forEach((d, c) => (ao[f + c] = radial * d * n[1 + c]));
    else [dx * dx, dx * dy, dx * dz, dy * dy, dy * dz, dz * dz].forEach((d, c) => (ao[f + c] = radial * d * n[4 + c]));
  }
  const mo = (row: number) => {
    let s = 0;
    for (let i = 0; i < v.nbf; i++) s += ao[i] * v.coefs[row * v.nbf + i];
    return s;
  };
  let rho = 0;
  for (let k = 0; k < v.nocc; k++) rho += v.weights[k] * mo(k) ** 2;
  return { rho, psi: mo(v.nocc) };
}
