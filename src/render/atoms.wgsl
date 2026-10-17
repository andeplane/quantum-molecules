// Ball-and-stick atoms: spheres are ray-cast impostors, bonds are instanced cylinders.

struct Camera {
  view: mat4x4f,
  proj: mat4x4f,
}

@group(0) @binding(0) var<uniform> cam: Camera;

struct Sphere {
  @location(0) center: vec3f,
  @location(1) radius: f32,
  @location(2) color: vec3f,
  @location(3) highlight: f32,
}

struct SphereOut {
  @builtin(position) pos: vec4f,
  @location(0) viewPos: vec3f,
  @location(1) center: vec3f,
  @location(2) radius: f32,
  @location(3) color: vec3f,
  @location(4) highlight: f32,
}

@vertex
fn sphereVs(@builtin(vertex_index) vi: u32, s: Sphere) -> SphereOut {
  var corners = array<vec2f, 4>(vec2f(-1.0, -1.0), vec2f(1.0, -1.0), vec2f(-1.0, 1.0), vec2f(1.0, 1.0));
  let c = (cam.view * vec4f(s.center, 1.0)).xyz;
  // Grow the quad and push it towards the camera so the whole sphere is covered.
  let grow = 1.6;
  let vp = c + vec3f(corners[vi] * s.radius * grow, s.radius);
  var o: SphereOut;
  o.pos = cam.proj * vec4f(vp, 1.0);
  o.viewPos = vp;
  o.center = c;
  o.radius = s.radius;
  o.color = s.color;
  o.highlight = s.highlight;
  return o;
}

fn shade(n: vec3f, v: vec3f, base: vec3f) -> vec3f {
  let l = normalize(vec3f(0.4, 0.6, 1.0));
  let diff = max(dot(n, l), 0.0);
  let spec = pow(max(dot(n, normalize(l + v)), 0.0), 48.0);
  let rim = pow(1.0 - max(dot(n, v), 0.0), 3.0);
  return base * (0.22 + 0.78 * diff) + vec3f(0.35 * spec) + 0.15 * rim * base;
}

struct FragOut {
  @location(0) color: vec4f,
  @builtin(frag_depth) depth: f32,
}

@fragment
fn sphereFs(in: SphereOut) -> FragOut {
  let dir = normalize(in.viewPos);
  let b = dot(dir, in.center);
  let disc = b * b - dot(in.center, in.center) + in.radius * in.radius;
  if (disc < 0.0) {
    discard;
  }
  let t = b - sqrt(disc);
  let hit = dir * t;
  let n = normalize(hit - in.center);
  let clip = cam.proj * vec4f(hit, 1.0);
  var base = in.color;
  if (in.highlight > 0.5) {
    base = mix(base, vec3f(1.0, 0.85, 0.3), 0.45);
  }
  var o: FragOut;
  o.color = vec4f(shade(n, -dir, base), 1.0);
  o.depth = clip.z / clip.w;
  return o;
}

struct Bond {
  @location(0) a: vec3f,
  @location(1) radius: f32,
  @location(2) b: vec3f,
  @location(3) colorA: vec3f,
  @location(4) colorB: vec3f,
}

struct BondOut {
  @builtin(position) pos: vec4f,
  @location(0) normal: vec3f,
  @location(1) viewPos: vec3f,
  @location(2) along: f32,
  @location(3) colorA: vec3f,
  @location(4) colorB: vec3f,
}

// `ring` is (cos, sin, along) on a unit cylinder from along=0 to along=1.
@vertex
fn bondVs(@location(5) ring: vec3f, bond: Bond) -> BondOut {
  let axis = bond.b - bond.a;
  let w = normalize(axis);
  let helper = select(vec3f(1.0, 0.0, 0.0), vec3f(0.0, 1.0, 0.0), abs(w.x) > 0.9);
  let u = normalize(cross(w, helper));
  let v = cross(w, u);
  let n = u * ring.x + v * ring.y;
  let p = bond.a + axis * ring.z + n * bond.radius;
  var o: BondOut;
  let vp = cam.view * vec4f(p, 1.0);
  o.pos = cam.proj * vp;
  o.viewPos = vp.xyz;
  o.normal = (cam.view * vec4f(n, 0.0)).xyz;
  o.along = ring.z;
  o.colorA = bond.colorA;
  o.colorB = bond.colorB;
  return o;
}

@fragment
fn bondFs(in: BondOut) -> @location(0) vec4f {
  let base = select(in.colorB, in.colorA, in.along < 0.5);
  return vec4f(shade(normalize(in.normal), normalize(-in.viewPos), base * 0.9), 1.0);
}
