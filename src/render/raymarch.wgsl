// Ray-marches the density cloud and the orbital isosurfaces, stopping at the atoms' depth.

struct Uniforms {
  invViewProj: mat4x4f,
  eye: vec3f,
  cloudOpacity: f32,
  boxMin: vec3f,
  logMin: f32,
  boxMax: vec3f,
  logMax: f32,
  voxel: vec3f,
  iso: f32,
  showCloud: f32,
  showOrbital: f32,
  orbitalOpacity: f32,
  frame: f32,
}

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var volume: texture_3d<f32>;
@group(0) @binding(2) var samp: sampler;
@group(0) @binding(3) var depthTex: texture_depth_2d;

struct VsOut {
  @builtin(position) pos: vec4f,
  @location(0) ndc: vec2f,
}

@vertex
fn vs(@builtin(vertex_index) i: u32) -> VsOut {
  let xy = vec2f(f32((i << 1u) & 2u), f32(i & 2u)) * 2.0 - 1.0;
  var o: VsOut;
  o.pos = vec4f(xy, 0.0, 1.0);
  o.ndc = xy;
  return o;
}

fn unproject(ndc: vec2f, z: f32) -> vec3f {
  let p = u.invViewProj * vec4f(ndc, z, 1.0);
  return p.xyz / p.w;
}

fn toUvw(p: vec3f) -> vec3f {
  return (p - u.boxMin) / (u.boxMax - u.boxMin);
}

fn sampleVol(p: vec3f) -> vec4f {
  return textureSampleLevel(volume, samp, toUvw(p), 0.0);
}

// A muted violet-to-white ramp for the cloud, kept apart from the blue and orange orbital lobes.
// Keep in sync with CLOUD_RAMP in src/ui/legend.ts.
fn cloudColor(t: f32) -> vec3f {
  let c0 = vec3f(0.30, 0.28, 0.40);
  let c1 = vec3f(0.66, 0.64, 0.78);
  let c2 = vec3f(1.00, 0.98, 0.95);
  return select(mix(c1, c2, (t - 0.5) * 2.0), mix(c0, c1, t * 2.0), t < 0.5);
}

fn hash(p: vec2f) -> f32 {
  return fract(sin(dot(p, vec2f(12.9898, 78.233)) + u.frame * 0.618) * 43758.5453);
}

fn orbitalNormal(p: vec3f) -> vec3f {
  let h = u.voxel;
  let g = vec3f(
    sampleVol(p + vec3f(h.x, 0.0, 0.0)).g - sampleVol(p - vec3f(h.x, 0.0, 0.0)).g,
    sampleVol(p + vec3f(0.0, h.y, 0.0)).g - sampleVol(p - vec3f(0.0, h.y, 0.0)).g,
    sampleVol(p + vec3f(0.0, 0.0, h.z)).g - sampleVol(p - vec3f(0.0, 0.0, h.z)).g,
  ) / (2.0 * h);
  return normalize(g + vec3f(1e-12));
}

@fragment
fn fs(in: VsOut) -> @location(0) vec4f {
  let frag = in.pos;
  let far = unproject(in.ndc, 1.0);
  let dir = normalize(far - u.eye);
  var tmax = 1e9;
  let depth = textureLoad(depthTex, vec2i(frag.xy), 0);
  if (depth < 1.0) {
    tmax = distance(unproject(in.ndc, depth), u.eye);
  }
  let inv = 1.0 / dir;
  let ta = (u.boxMin - u.eye) * inv;
  let tb = (u.boxMax - u.eye) * inv;
  let t0 = max(max(max(min(ta.x, tb.x), min(ta.y, tb.y)), min(ta.z, tb.z)), 0.0);
  let t1 = min(min(min(max(ta.x, tb.x), max(ta.y, tb.y)), max(ta.z, tb.z)), tmax);
  if (t1 <= t0) {
    discard;
  }
  let dt = 0.6 * min(min(u.voxel.x, u.voxel.y), u.voxel.z);
  let steps = min(u32((t1 - t0) / dt) + 1u, 1024u);
  var t = t0 + dt * hash(frag.xy);
  var color = vec3f(0.0);
  var alpha = 0.0;
  var prevPsi = sampleVol(u.eye + dir * t).g;
  let light = normalize(-dir + vec3f(0.3, 0.5, 0.0));
  for (var i = 0u; i < steps; i++) {
    if (t > t1 || alpha > 0.985) {
      break;
    }
    let p = u.eye + dir * t;
    let s = sampleVol(p);
    if (u.showOrbital > 0.5) {
      let psi = s.g;
      // A crossing of +iso or -iso between the previous sample and this one.
      let crossPos = (prevPsi - u.iso) * (psi - u.iso) < 0.0;
      let crossNeg = (prevPsi + u.iso) * (psi + u.iso) < 0.0;
      if (crossPos || crossNeg) {
        let level = select(-u.iso, u.iso, crossPos);
        let f = clamp((level - prevPsi) / (psi - prevPsi), 0.0, 1.0);
        let hit = p - dir * dt * (1.0 - f);
        var n = orbitalNormal(hit);
        if (dot(n, dir) > 0.0) {
          n = -n;
        }
        // Keep in sync with LOBE_POSITIVE and LOBE_NEGATIVE in src/ui/legend.ts.
        let base = select(vec3f(0.95, 0.42, 0.18), vec3f(0.22, 0.52, 1.0), level > 0.0);
        let diff = max(dot(n, light), 0.0);
        let spec = pow(max(dot(n, normalize(light - dir)), 0.0), 40.0);
        let shaded = base * (0.25 + 0.75 * diff) + vec3f(0.35 * spec);
        let a = u.orbitalOpacity;
        color += (1.0 - alpha) * a * shaded;
        alpha += (1.0 - alpha) * a;
      }
      prevPsi = psi;
    }
    if (u.showCloud > 0.5) {
      let rho = max(s.r, 1e-12);
      let x = clamp((log2(rho) * 0.30103 - u.logMin) / (u.logMax - u.logMin), 0.0, 1.0);
      let a = 1.0 - exp(-u.cloudOpacity * x * x * x * dt);
      color += (1.0 - alpha) * a * cloudColor(x) * (0.8 + 0.8 * x);
      alpha += (1.0 - alpha) * a;
    }
    t += dt;
  }
  return vec4f(color, alpha);
}
