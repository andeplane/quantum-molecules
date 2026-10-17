// Evaluates the electron density and one molecular orbital on a 3D grid.
// The basis layout matches `gpu_basis` in crates/qchem-wasm and `evalAo` in src/chem/orbitals.ts.

const MAX_BF: u32 = 128u;

struct Params {
  boxMin: vec3f,
  nbf: u32,
  boxMax: vec3f,
  nshell: u32,
  dims: vec3u,
  // Rows of `coefs`: `nocc` occupied orbitals weighted by `weights`, then the selected orbital.
  nocc: u32,
}

struct Shell {
  center: vec3f,
  l: u32,
  nprim: u32,
  prim: u32,
  fn0: u32,
  // Smallest exponent, for skipping shells that vanish at this point.
  minExp: f32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> shells: array<Shell>;
@group(0) @binding(2) var<storage, read> prims: array<vec2f>;
// Cartesian component norms: [s, px, py, pz, dxx, dxy, dxz, dyy, dyz, dzz].
@group(0) @binding(3) var<storage, read> norms: array<f32, 10>;
@group(0) @binding(4) var<storage, read> coefs: array<f32>;
@group(0) @binding(5) var<storage, read> weights: array<f32>;
@group(0) @binding(6) var volume: texture_storage_3d<rgba16float, write>;

var<private> ao: array<f32, MAX_BF>;

fn evalAo(p: vec3f) {
  for (var s = 0u; s < params.nshell; s++) {
    let sh = shells[s];
    let d = p - sh.center;
    let r2 = dot(d, d);
    let f = sh.fn0;
    var radial = 0.0;
    if (sh.minExp * r2 < 40.0) {
      for (var k = 0u; k < sh.nprim; k++) {
        let pr = prims[sh.prim + k];
        radial += pr.y * exp(-pr.x * r2);
      }
    }
    if (sh.l == 0u) {
      ao[f] = radial * norms[0];
    } else if (sh.l == 1u) {
      ao[f] = radial * d.x * norms[1];
      ao[f + 1u] = radial * d.y * norms[2];
      ao[f + 2u] = radial * d.z * norms[3];
    } else {
      ao[f] = radial * d.x * d.x * norms[4];
      ao[f + 1u] = radial * d.x * d.y * norms[5];
      ao[f + 2u] = radial * d.x * d.z * norms[6];
      ao[f + 3u] = radial * d.y * d.y * norms[7];
      ao[f + 4u] = radial * d.y * d.z * norms[8];
      ao[f + 5u] = radial * d.z * d.z * norms[9];
    }
  }
}

fn moValue(row: u32) -> f32 {
  var v = 0.0;
  let base = row * params.nbf;
  for (var i = 0u; i < params.nbf; i++) {
    v += ao[i] * coefs[base + i];
  }
  return v;
}

@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) id: vec3u) {
  if (any(id >= params.dims)) {
    return;
  }
  let t = vec3f(id) / vec3f(max(params.dims - 1u, vec3u(1u)));
  let p = mix(params.boxMin, params.boxMax, t);
  evalAo(p);
  var rho = 0.0;
  for (var k = 0u; k < params.nocc; k++) {
    let v = moValue(k);
    rho += weights[k] * v * v;
  }
  let psi = moValue(params.nocc);
  textureStore(volume, id, vec4f(rho, psi, 0.0, 1.0));
}
