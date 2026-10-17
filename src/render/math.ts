/** Small column-major 4x4 matrix and vector helpers (the layout WGSL expects). */

export type Vec3 = [number, number, number];
export type Mat4 = Float32Array<ArrayBuffer>;

export const sub = (a: Vec3, b: Vec3): Vec3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const scale = (a: Vec3, s: number): Vec3 => [a[0] * s, a[1] * s, a[2] * s];
export const dot = (a: Vec3, b: Vec3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a: Vec3, b: Vec3): Vec3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const length = (a: Vec3) => Math.hypot(a[0], a[1], a[2]);
export const normalize = (a: Vec3): Vec3 => scale(a, 1 / (length(a) || 1));

export function multiply(a: Mat4, b: Mat4): Mat4 {
  const out = new Float32Array(16);
  for (let c = 0; c < 4; c++)
    for (let r = 0; r < 4; r++) {
      let s = 0;
      for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
      out[c * 4 + r] = s;
    }
  return out;
}

/** Perspective projection with WebGPU's [0, 1] depth range. */
export function perspective(fovy: number, aspect: number, near: number, far: number): Mat4 {
  const f = 1 / Math.tan(fovy / 2);
  const out = new Float32Array(16);
  out[0] = f / aspect;
  out[5] = f;
  out[10] = far / (near - far);
  out[11] = -1;
  out[14] = (far * near) / (near - far);
  return out;
}

export function lookAt(eye: Vec3, target: Vec3, up: Vec3): Mat4 {
  const z = normalize(sub(eye, target));
  const x = normalize(cross(up, z));
  const y = cross(z, x);
  // prettier-ignore
  return new Float32Array([
    x[0], y[0], z[0], 0,
    x[1], y[1], z[1], 0,
    x[2], y[2], z[2], 0,
    -dot(x, eye), -dot(y, eye), -dot(z, eye), 1,
  ]);
}

export function invert(m: Mat4): Mat4 {
  const inv = new Float64Array(16);
  const a = Array.from(m);
  // Gauss-Jordan on [m | I], rows of the column-major matrix accessed as a[c*4+r].
  const w = Array.from({ length: 4 }, (_, r) => [0, 1, 2, 3].map((c) => a[c * 4 + r]).concat([0, 1, 2, 3].map((c) => (c === r ? 1 : 0))));
  for (let col = 0; col < 4; col++) {
    let piv = col;
    for (let r = col + 1; r < 4; r++) if (Math.abs(w[r][col]) > Math.abs(w[piv][col])) piv = r;
    [w[col], w[piv]] = [w[piv], w[col]];
    const p = w[col][col];
    if (Math.abs(p) < 1e-30) throw new Error('singular matrix');
    for (let c = 0; c < 8; c++) w[col][c] /= p;
    for (let r = 0; r < 4; r++) {
      if (r === col) continue;
      const f = w[r][col];
      for (let c = 0; c < 8; c++) w[r][c] -= f * w[col][c];
    }
  }
  for (let r = 0; r < 4; r++) for (let c = 0; c < 4; c++) inv[c * 4 + r] = w[r][4 + c];
  return new Float32Array(inv);
}

export function transformPoint(m: Mat4, p: Vec3): Vec3 {
  const x = m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12];
  const y = m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13];
  const z = m[2] * p[0] + m[6] * p[1] + m[10] * p[2] + m[14];
  const w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
  return [x / w, y / w, z / w];
}

/** Distance along a ray to a sphere, or null if it misses. */
export function raySphere(origin: Vec3, dir: Vec3, center: Vec3, radius: number): number | null {
  const oc = sub(origin, center);
  const b = dot(oc, dir);
  const c = dot(oc, oc) - radius * radius;
  const disc = b * b - c;
  if (disc < 0) return null;
  const t = -b - Math.sqrt(disc);
  return t > 0 ? t : null;
}
