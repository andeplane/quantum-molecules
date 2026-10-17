(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const n of document.querySelectorAll('link[rel="modulepreload"]'))o(n);new MutationObserver(n=>{for(const i of n)if(i.type==="childList")for(const l of i.addedNodes)l.tagName==="LINK"&&l.rel==="modulepreload"&&o(l)}).observe(document,{childList:!0,subtree:!0});function r(n){const i={};return n.integrity&&(i.integrity=n.integrity),n.referrerPolicy&&(i.referrerPolicy=n.referrerPolicy),n.crossOrigin==="use-credentials"?i.credentials="include":n.crossOrigin==="anonymous"?i.credentials="omit":i.credentials="same-origin",i}function o(n){if(n.ep)return;n.ep=!0;const i=r(n);fetch(n.href,i)}})();const D=1/.52917721092,k=27.211386245988,Se=2.541746473,Ae=[["H",1.008,.31,[.95,.95,.95]],["He",4.0026,.28,[.85,1,1]],["Li",6.94,1.28,[.8,.5,1]],["Be",9.0122,.96,[.76,1,0]],["B",10.81,.84,[1,.71,.71]],["C",12.011,.76,[.38,.38,.38]],["N",14.007,.71,[.19,.31,.97]],["O",15.999,.66,[1,.13,.13]],["F",18.998,.57,[.56,.88,.31]],["Ne",20.18,.58,[.7,.89,.96]]],Ce=Ae.map(([t,e,r,o],n)=>({symbol:t,z:n+1,mass:e,covalentRadius:r,color:o})),Oe=new Map(Ce.map(t=>[t.symbol.toLowerCase(),t]));function N(t){const e=Oe.get(t.toLowerCase());if(!e)throw new Error(`Unsupported element: ${t} (only H to Ne are available)`);return e}function Ue(t){return N(t.symbol).z}function ce(t){return t.atoms.reduce((e,r)=>e+Ue(r),0)-t.charge}function Le(t,e=1.25){const r=[];for(let o=0;o<t.length;o++)for(let n=0;n<o;n++){const[i,l,c]=t[o].xyz,[f,u,d]=t[n].xyz;Math.hypot(i-f,l-u,c-d)/D<e*(N(t[o].symbol).covalentRadius+N(t[n].symbol).covalentRadius)&&r.push([n,o])}return r}function ue(t){const e=[];for(const r of t.split(/\r?\n/)){const o=r.trim().split(/\s+/);if(o.length<4)continue;const[n,...i]=o,l=i.slice(0,3).map(Number);if(!/^[A-Za-z]{1,2}$/.test(n)||l.some(f=>!Number.isFinite(f)))continue;const c=n[0].toUpperCase()+n.slice(1).toLowerCase();N(c),e.push({symbol:c,xyz:l.map(f=>f*D)})}return e}function Fe(t,e=""){const r=t.map(o=>`${o.symbol.padEnd(2)} ${o.xyz.map(n=>(n/D).toFixed(6).padStart(12)).join(" ")}`);return`${t.length}
${e}
${r.join(`
`)}
`}const E=["sto-3g","6-31g","6-31g*"];function Re(){const t=[];for(let e=0;e<6;e++){const r=e*Math.PI/3;t.push(`C ${(1.397*Math.cos(r)).toFixed(4)} ${(1.397*Math.sin(r)).toFixed(4)} 0`)}for(let e=0;e<6;e++){const r=e*Math.PI/3;t.push(`H ${(2.481*Math.cos(r)).toFixed(4)} ${(2.481*Math.sin(r)).toFixed(4)} 0`)}return t.join(`
`)}const fe=[{id:"h2",name:"Hydrogen H₂",xyz:`H 0 0 -0.3705
H 0 0 0.3705`,charge:0,multiplicity:1,bases:E},{id:"h2o",name:"Water H₂O",xyz:`O 0 0 0.1173
H 0 0.7572 -0.4692
H 0 -0.7572 -0.4692`,charge:0,multiplicity:1,bases:E},{id:"nh3",name:"Ammonia NH₃",xyz:`N 0 0 0.1162
H 0 0.9377 -0.2711
H 0.8121 -0.4689 -0.2711
H -0.8121 -0.4689 -0.2711`,charge:0,multiplicity:1,bases:E},{id:"ch4",name:"Methane CH₄",xyz:`C 0 0 0
H 0.6291 0.6291 0.6291
H -0.6291 -0.6291 0.6291
H -0.6291 0.6291 -0.6291
H 0.6291 -0.6291 -0.6291`,charge:0,multiplicity:1,bases:E},{id:"hf",name:"Hydrogen fluoride HF",xyz:`F 0 0 0.0917
H 0 0 -0.8253`,charge:0,multiplicity:1,bases:E},{id:"n2",name:"Nitrogen N₂",xyz:`N 0 0 -0.5488
N 0 0 0.5488`,charge:0,multiplicity:1,bases:E},{id:"co2",name:"Carbon dioxide CO₂",xyz:`C 0 0 0
O 0 0 1.1621
O 0 0 -1.1621`,charge:0,multiplicity:1,bases:E},{id:"ch2o",name:"Formaldehyde CH₂O",xyz:`C 0 0 -0.5285
O 0 0 0.6745
H 0 0.9363 -1.1079
H 0 -0.9363 -1.1079`,charge:0,multiplicity:1,bases:E},{id:"c2h4",name:"Ethylene C₂H₄",xyz:`C 0 0 0.6695
C 0 0 -0.6695
H 0 0.9288 1.2321
H 0 -0.9288 1.2321
H 0 0.9288 -1.2321
H 0 -0.9288 -1.2321`,charge:0,multiplicity:1,bases:E},{id:"c6h6",name:"Benzene C₆H₆",xyz:Re(),charge:0,multiplicity:1,bases:["sto-3g"]},{id:"oh",name:"Hydroxyl radical OH",xyz:`O 0 0 0
H 0 0 0.9697`,charge:0,multiplicity:2,bases:E},{id:"o2",name:"Oxygen O₂ (triplet)",xyz:`O 0 0 -0.6038
O 0 0 0.6038`,charge:0,multiplicity:3,bases:E}];function Ge(t){const e=fe.find(r=>r.id===t);if(!e)throw new Error(`Unknown preset ${t}`);return e}function de(t,e="alpha"){const r=e==="beta"&&t.unrestricted?t.epsB:t.epsA,o=e==="beta"?t.nbeta:t.nalpha;return r.map((n,i)=>{let l="";return i===o-1?l="HOMO":i===o?l="LUMO":i<o?l=`HOMO−${o-1-i}`:l=`LUMO+${i-o}`,{index:i,energy:n,occupied:i<o,label:l}})}function Ne(t){let e=t.epsA[t.nalpha-1]??-1/0,r=t.epsA[t.nalpha]??1/0;return t.unrestricted&&(t.nbeta>0&&(e=Math.max(e,t.epsB[t.nbeta-1])),r=Math.min(r,t.epsB[t.nbeta]??1/0)),{homo:e,lumo:Number.isFinite(r)?r:null}}class Te{worker;nextId=1;pending=new Map;relaxes=new Map;busy=!1;queued=null;ready;constructor(e){this.worker=e;let r;this.ready=new Promise(o=>r=o),e.onmessage=o=>{const n=o.data;if(n.type==="ready")return r();if(n.type==="relax-step"||n.type==="relax-done"){const l=this.relaxes.get(n.id);if(!l)return;n.type==="relax-step"?l.onStep(n.output,n.step):(this.relaxes.delete(n.id),l.onDone(n.converged,n.steps));return}if(n.type==="error"&&this.relaxes.has(n.id)){const l=this.relaxes.get(n.id);this.relaxes.delete(n.id),l.onError(n.message);return}const i=this.pending.get(n.id);i&&(this.pending.delete(n.id),n.type==="result"?i.resolve(n.output):n.type==="error"&&i.reject(new Error(n.message)))}}send(e){this.worker.postMessage(e)}compute(e){const r=this.nextId++;return new Promise((o,n)=>{this.pending.set(r,{resolve:o,reject:n}),this.send({id:r,type:"compute",input:e})})}computeLatest(e){return this.busy?(this.queued?.resolve(null),new Promise(r=>this.queued={input:e,resolve:r})):(this.busy=!0,this.compute(e).finally(()=>{this.busy=!1;const r=this.queued;this.queued=null,r&&this.computeLatest(r.input).then(r.resolve,()=>r.resolve(null))}))}relax(e,r,o=60){const n=this.nextId++;return this.relaxes.set(n,r),this.send({id:n,type:"relax",input:e,maxSteps:o}),n}cancelRelax(e){this.relaxes.delete(e),this.send({id:e,type:"cancel"})}reset(){this.send({id:this.nextId++,type:"reset"})}}const q=(t,e)=>[t[0]-e[0],t[1]-e[1],t[2]-e[2]],he=(t,e)=>[t[0]+e[0],t[1]+e[1],t[2]+e[2]],K=(t,e)=>[t[0]*e,t[1]*e,t[2]*e],L=(t,e)=>t[0]*e[0]+t[1]*e[1]+t[2]*e[2],oe=(t,e)=>[t[1]*e[2]-t[2]*e[1],t[2]*e[0]-t[0]*e[2],t[0]*e[1]-t[1]*e[0]],$e=t=>Math.hypot(t[0],t[1],t[2]),Y=t=>K(t,1/($e(t)||1));function He(t,e){const r=new Float32Array(16);for(let o=0;o<4;o++)for(let n=0;n<4;n++){let i=0;for(let l=0;l<4;l++)i+=t[l*4+n]*e[o*4+l];r[o*4+n]=i}return r}function Ie(t,e,r,o){const n=1/Math.tan(t/2),i=new Float32Array(16);return i[0]=n/e,i[5]=n,i[10]=o/(r-o),i[11]=-1,i[14]=o*r/(r-o),i}function _e(t,e,r){const o=Y(q(t,e)),n=Y(oe(r,o)),i=oe(o,n);return new Float32Array([n[0],i[0],o[0],0,n[1],i[1],o[1],0,n[2],i[2],o[2],0,-L(n,t),-L(i,t),-L(o,t),1])}function pe(t){const e=new Float64Array(16),r=Array.from(t),o=Array.from({length:4},(n,i)=>[0,1,2,3].map(l=>r[l*4+i]).concat([0,1,2,3].map(l=>l===i?1:0)));for(let n=0;n<4;n++){let i=n;for(let c=n+1;c<4;c++)Math.abs(o[c][n])>Math.abs(o[i][n])&&(i=c);[o[n],o[i]]=[o[i],o[n]];const l=o[n][n];if(Math.abs(l)<1e-30)throw new Error("singular matrix");for(let c=0;c<8;c++)o[n][c]/=l;for(let c=0;c<4;c++){if(c===n)continue;const f=o[c][n];for(let u=0;u<8;u++)o[c][u]-=f*o[n][u]}}for(let n=0;n<4;n++)for(let i=0;i<4;i++)e[i*4+n]=o[n][4+i];return new Float32Array(e)}function ie(t,e){const r=t[0]*e[0]+t[4]*e[1]+t[8]*e[2]+t[12],o=t[1]*e[0]+t[5]*e[1]+t[9]*e[2]+t[13],n=t[2]*e[0]+t[6]*e[1]+t[10]*e[2]+t[14],i=t[3]*e[0]+t[7]*e[1]+t[11]*e[2]+t[15];return[r/i,o/i,n/i]}function Ve(t,e,r,o){const n=q(t,r),i=L(n,e),l=L(n,n)-o*o,c=i*i-l;if(c<0)return null;const f=-i-Math.sqrt(c);return f>0?f:null}class ke{target=[0,0,0];distance=15;yaw=.6;pitch=.35;fovy=40*Math.PI/180;aspect=1;get eye(){const e=Math.cos(this.pitch);return he(this.target,K([e*Math.sin(this.yaw),Math.sin(this.pitch),e*Math.cos(this.yaw)],this.distance))}get view(){return _e(this.eye,this.target,[0,1,0])}get projection(){return Ie(this.fovy,this.aspect,Math.max(.05,this.distance*.02),this.distance*10+50)}get viewProjection(){return He(this.projection,this.view)}orbit(e,r){this.yaw-=e*.008,this.pitch=Math.max(-1.5,Math.min(1.5,this.pitch+r*.008))}zoom(e){this.distance=Math.max(2,Math.min(200,this.distance*e))}fit(e,r){this.target=e,this.distance=Math.max(6,(r+2.5)/Math.sin(this.fovy/2))}ray(e,r){const o=pe(this.viewProjection),n=ie(o,[e,r,0]),i=ie(o,[e,r,1]);return{origin:n,dir:Y(q(i,n))}}basis(){const e=this.view;return{right:[e[0],e[4],e[8]],up:[e[1],e[5],e[9]],forward:[-e[2],-e[6],-e[10]]}}}const De=`// Ball-and-stick atoms: spheres are ray-cast impostors, bonds are instanced cylinders.

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

// \`ring\` is (cos, sin, along) on a unit cylinder from along=0 to along=1.
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
`,qe=`// Ray-marches the density cloud and the orbital isosurfaces, stopping at the atoms' depth.

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

// A perceptual blue-to-white ramp for the cloud.
fn cloudColor(t: f32) -> vec3f {
  let c0 = vec3f(0.10, 0.22, 0.55);
  let c1 = vec3f(0.25, 0.65, 0.95);
  let c2 = vec3f(0.95, 0.98, 1.00);
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
`,je=`// Evaluates the electron density and one molecular orbital on a 3D grid.
// The basis layout matches \`gpu_basis\` in crates/qchem-wasm and \`evalAo\` in src/chem/orbitals.ts.

const MAX_BF: u32 = 128u;

struct Params {
  boxMin: vec3f,
  nbf: u32,
  boxMax: vec3f,
  nshell: u32,
  dims: vec3u,
  // Rows of \`coefs\`: \`nocc\` occupied orbitals weighted by \`weights\`, then the selected orbital.
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
`,ae="depth24plus",j=18;function We(){const t=[];for(let e=0;e<j;e++){const r=e/j*2*Math.PI,o=(e+1)/j*2*Math.PI,n=(i,l)=>t.push(Math.cos(i),Math.sin(i),l);n(r,0),n(o,0),n(o,1),n(r,0),n(o,1),n(r,1)}return new Float32Array(t)}class J{device;context;format;depth=null;cameraBuffer;rayBuffer;paramsBuffer;normsBuffer;ringBuffer;ringCount;sphereBuffer=null;sphereCount=0;bondBuffer=null;bondCount=0;spherePipeline;bondPipeline;rayPipeline;volumePipeline;atomsBindGroup;sampler;volume=null;box=null;frame=0;adapter=null;background={r:.035,g:.045,b:.07,a:1};static async create(e){if(!navigator.gpu)throw new Error("WebGPU is not available in this browser.");const r=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!r)throw new Error("No WebGPU adapter was found.");const o=await r.requestDevice(),n=e?e.getContext("webgpu"):null;if(e&&!n)throw new Error("Could not create a WebGPU canvas context.");const i=new J(o,n,navigator.gpu.getPreferredCanvasFormat());return i.adapter=r,i}constructor(e,r,o){this.device=e,this.context=r,this.format=o,r?.configure({device:e,format:o,alphaMode:"opaque"});const n=GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST;this.cameraBuffer=e.createBuffer({size:128,usage:n}),this.rayBuffer=e.createBuffer({size:144,usage:n}),this.paramsBuffer=e.createBuffer({size:48,usage:n}),this.normsBuffer=e.createBuffer({size:40,usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});const i=We();this.ringCount=i.length/3,this.ringBuffer=e.createBuffer({size:i.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),e.queue.writeBuffer(this.ringBuffer,0,i),this.sampler=e.createSampler({magFilter:"linear",minFilter:"linear"});const l=e.createShaderModule({code:De,label:"atoms"}),c=e.createPipelineLayout({bindGroupLayouts:[this.cameraLayout()]}),f={format:ae,depthWriteEnabled:!0,depthCompare:"less"},u=(m,h)=>({shaderLocation:m,offset:h,format:"float32x3"}),d=(m,h)=>({shaderLocation:m,offset:h,format:"float32"});this.spherePipeline=e.createRenderPipeline({label:"spheres",layout:c,vertex:{module:l,entryPoint:"sphereVs",buffers:[{arrayStride:32,stepMode:"instance",attributes:[u(0,0),d(1,12),u(2,16),d(3,28)]}]},fragment:{module:l,entryPoint:"sphereFs",targets:[{format:o}]},primitive:{topology:"triangle-strip"},depthStencil:f}),this.bondPipeline=e.createRenderPipeline({label:"bonds",layout:c,vertex:{module:l,entryPoint:"bondVs",buffers:[{arrayStride:52,stepMode:"instance",attributes:[u(0,0),d(1,12),u(2,16),u(3,28),u(4,40)]},{arrayStride:12,stepMode:"vertex",attributes:[u(5,0)]}]},fragment:{module:l,entryPoint:"bondFs",targets:[{format:o}]},primitive:{topology:"triangle-list",cullMode:"back"},depthStencil:f}),this.atomsBindGroup=e.createBindGroup({layout:this.cameraLayout(),entries:[{binding:0,resource:{buffer:this.cameraBuffer}}]});const p=e.createShaderModule({code:qe,label:"raymarch"});this.rayPipeline=e.createRenderPipeline({label:"raymarch",layout:"auto",vertex:{module:p,entryPoint:"vs"},fragment:{module:p,entryPoint:"fs",targets:[{format:o,blend:{color:{srcFactor:"one",dstFactor:"one-minus-src-alpha"},alpha:{srcFactor:"one",dstFactor:"one-minus-src-alpha"}}}]}}),this.volumePipeline=e.createComputePipeline({label:"volume",layout:"auto",compute:{module:e.createShaderModule({code:je,label:"volume"}),entryPoint:"main"}})}cameraLayoutCache=null;cameraLayout(){return this.cameraLayoutCache??=this.device.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.VERTEX|GPUShaderStage.FRAGMENT,buffer:{type:"uniform"}}]}),this.cameraLayoutCache}resize(e,r){this.depth&&this.depth.width===e&&this.depth.height===r||(this.depth?.destroy(),this.depth=this.device.createTexture({size:[e,r],format:ae,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.TEXTURE_BINDING}))}setAtoms(e,r){const o=new Float32Array(Math.max(1,e.length)*8);e.forEach((i,l)=>o.set([...i.center,i.radius,...i.color,i.highlight?1:0],8*l)),(!this.sphereBuffer||this.sphereBuffer.size<o.byteLength)&&(this.sphereBuffer?.destroy(),this.sphereBuffer=this.device.createBuffer({size:o.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST})),this.device.queue.writeBuffer(this.sphereBuffer,0,o),this.sphereCount=e.length;const n=new Float32Array(Math.max(1,r.length)*13);r.forEach((i,l)=>n.set([...i.a,i.radius,...i.b,...i.colorA,...i.colorB],13*l)),(!this.bondBuffer||this.bondBuffer.size<n.byteLength)&&(this.bondBuffer?.destroy(),this.bondBuffer=this.device.createBuffer({size:n.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST})),this.device.queue.writeBuffer(this.bondBuffer,0,n),this.bondCount=r.length}updateVolume(e){const{device:r}=this,[o,n,i]=e.box.dims;(!this.volume||this.volume.width!==o||this.volume.height!==n||this.volume.depthOrArrayLayers!==i)&&(this.volume?.destroy(),this.volume=r.createTexture({size:[o,n,i],dimension:"3d",format:"rgba16float",usage:GPUTextureUsage.STORAGE_BINDING|GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_SRC})),this.box=e.box;const l=new ArrayBuffer(48),c=new Float32Array(l),f=new Uint32Array(l);c.set(e.box.min,0),f[3]=e.nbf,c.set(e.box.max,4),f[7]=e.nshell,f.set(e.box.dims,8),f[11]=e.nocc,r.queue.writeBuffer(this.paramsBuffer,0,l),r.queue.writeBuffer(this.normsBuffer,0,e.norms);const u=b=>{const x=r.createBuffer({size:Math.max(16,b.byteLength),usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});return r.queue.writeBuffer(x,0,b),x},d=[u(e.shells),u(e.prims),u(e.coefs),u(e.weights.length?e.weights:new Float32Array(4))],p=r.createBindGroup({layout:this.volumePipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.paramsBuffer}},{binding:1,resource:{buffer:d[0]}},{binding:2,resource:{buffer:d[1]}},{binding:3,resource:{buffer:this.normsBuffer}},{binding:4,resource:{buffer:d[2]}},{binding:5,resource:{buffer:d[3]}},{binding:6,resource:this.volume.createView()}]}),m=r.createCommandEncoder(),h=m.beginComputePass();h.setPipeline(this.volumePipeline),h.setBindGroup(0,p),h.dispatchWorkgroups(Math.ceil(o/4),Math.ceil(n/4),Math.ceil(i/4)),h.end(),r.queue.submit([m.finish()]);for(const b of d)b.destroy()}render(e,r,o){const{device:n}=this,i=o??this.context?.getCurrentTexture();if(!i)return;this.resize(i.width,i.height),e.aspect=i.width/i.height;const l=new Float32Array(32);l.set(e.view,0),l.set(e.projection,16),n.queue.writeBuffer(this.cameraBuffer,0,l);const c=n.createCommandEncoder(),f=c.beginRenderPass({colorAttachments:[{view:i.createView(),clearValue:this.background,loadOp:"clear",storeOp:"store"}],depthStencilAttachment:{view:this.depth.createView(),depthClearValue:1,depthLoadOp:"clear",depthStoreOp:"store"}});if(this.bondCount>0&&this.bondBuffer&&(f.setPipeline(this.bondPipeline),f.setBindGroup(0,this.atomsBindGroup),f.setVertexBuffer(0,this.bondBuffer),f.setVertexBuffer(1,this.ringBuffer),f.draw(this.ringCount,this.bondCount)),this.sphereCount>0&&this.sphereBuffer&&(f.setPipeline(this.spherePipeline),f.setBindGroup(0,this.atomsBindGroup),f.setVertexBuffer(0,this.sphereBuffer),f.draw(4,this.sphereCount)),f.end(),this.volume&&this.box&&(r.showCloud||r.showOrbital)){const u=this.box,d=new Float32Array(36);d.set(pe(e.viewProjection),0),d.set(e.eye,16),d[19]=r.cloudOpacity,d.set(u.min,20),d[23]=r.logMin,d.set(u.max,24),d[27]=r.logMax,d.set(u.max.map((h,b)=>(h-u.min[b])/Math.max(1,u.dims[b]-1)),28),d[31]=r.iso,d[32]=r.showCloud?1:0,d[33]=r.showOrbital?1:0,d[34]=r.orbitalOpacity,d[35]=this.frame++%64,n.queue.writeBuffer(this.rayBuffer,0,d);const p=n.createBindGroup({layout:this.rayPipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.rayBuffer}},{binding:1,resource:this.volume.createView()},{binding:2,resource:this.sampler},{binding:3,resource:this.depth.createView()}]}),m=c.beginRenderPass({colorAttachments:[{view:i.createView(),loadOp:"load",storeOp:"store"}]});m.setPipeline(this.rayPipeline),m.setBindGroup(0,p),m.draw(3),m.end()}n.queue.submit([c.finish()])}async snapshot(e,r,o,n){const i=this.device.createTexture({size:[o,n],format:this.format,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});this.render(e,r,i);const l=Math.ceil(o*4/256)*256,c=this.device.createBuffer({size:l*n,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),f=this.device.createCommandEncoder();f.copyTextureToBuffer({texture:i},{buffer:c,bytesPerRow:l},[o,n]),this.device.queue.submit([f.finish()]),await c.mapAsync(GPUMapMode.READ);const u=new Uint8Array(c.getMappedRange()),d=new Uint8ClampedArray(o*n*4),p=this.format.startsWith("bgra");for(let m=0;m<n;m++)for(let h=0;h<o;h++){const b=m*l+4*h,x=4*(m*o+h);d[x]=u[b+(p?2:0)],d[x+1]=u[b+1],d[x+2]=u[b+(p?0:2)],d[x+3]=255}return c.unmap(),c.destroy(),i.destroy(),d}async readVolume(){if(!this.volume||!this.box)return null;const[e,r,o]=this.box.dims,n=Math.ceil(e*8/256)*256,i=this.device.createBuffer({size:n*r*o,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),l=this.device.createCommandEncoder();l.copyTextureToBuffer({texture:this.volume},{buffer:i,bytesPerRow:n,rowsPerImage:r},[e,r,o]),this.device.queue.submit([l.finish()]),await i.mapAsync(GPUMapMode.READ);const c=new Uint16Array(i.getMappedRange()),f=new Float32Array(e*r*o*2);for(let u=0;u<o;u++)for(let d=0;d<r;d++)for(let p=0;p<e;p++){const m=(u*r+d)*n/2+p*4,h=2*((u*r+d)*e+p);f[h]=se(c[m]),f[h+1]=se(c[m+1])}return i.unmap(),i.destroy(),{dims:[e,r,o],data:f}}}function se(t){const e=t&32768?-1:1,r=t>>10&31,o=t&1023;return r===0?e*o*2**-24:r===31?o?NaN:e*(1/0):e*(1+o/1024)*2**(r-15)}const le=128,Xe=5;function Ye(t,e=128,r=Xe){const o=[1/0,1/0,1/0],n=[-1/0,-1/0,-1/0];for(let f=0;f<t.length;f+=3)for(let u=0;u<3;u++)o[u]=Math.min(o[u],t[f+u]-r),n[u]=Math.max(n[u],t[f+u]+r);const i=n.map((f,u)=>f-o[u]),l=Math.max(...i)/(e-1),c=i.map(f=>Math.min(e,Math.ceil(f/l-1e-9)+1));for(let f=0;f<3;f++){const u=(o[f]+n[f])/2,d=(c[f]-1)*l/2;o[f]=u-d,n[f]=u+d}return{min:o,max:n,dims:c}}function Ze(t,e,r=t,o=128){const{nbf:n,nmo:i}=t;if(n>le)throw new Error(`The renderer supports at most ${le} basis functions (this system has ${n}).`);const{shells:l,prims:c,componentNorms:f}=t.basis,u=l.length/8,d=new ArrayBuffer(32*u),p=new Float32Array(d),m=new Uint32Array(d);for(let v=0;v<u;v++){const g=8*v,z=l[g+4],w=l[g+5];let G=1/0;for(let B=0;B<z;B++)G=Math.min(G,c[2*(w+B)]);p[g]=l[g],p[g+1]=l[g+1],p[g+2]=l[g+2],m[g+3]=l[g+3],m[g+4]=z,m[g+5]=w,m[g+6]=l[g+6],p[g+7]=G}const h=t.unrestricted?r.cb:r.ca,b=[];for(let v=0;v<t.nalpha;v++)b.push({c:r.ca,k:v,w:t.unrestricted?1:2});if(t.unrestricted)for(let v=0;v<t.nbeta;v++)b.push({c:h,k:v,w:1});const x=b.length,C=new Float32Array((x+1)*n);if(b.forEach(({c:v,k:g},z)=>{for(let w=0;w<n;w++)C[z*n+w]=v[w*i+g]}),e){const v=e.spin==="beta"?h:r.ca;for(let g=0;g<n;g++)C[x*n+g]=v[g*i+e.index]}return{nbf:n,nshell:u,nocc:x,shells:d,prims:Float32Array.from(c),norms:Float32Array.from(f.flat()),coefs:C,weights:Float32Array.from(b.map(v=>v.w)),box:Ye(t.coordinates,o)}}const W=t=>`${(t*k).toFixed(2)} eV`;function Ke(t){const{homo:e,lumo:r}=Ne(t),o=Math.hypot(...t.dipole)*Se,n=[{label:"Total energy",value:`${t.energy.toFixed(6)} Ha`,title:`${(t.energy*k).toFixed(3)} eV`},{label:"HOMO",value:W(e)},{label:"LUMO",value:r===null?"none":W(r)},{label:"HOMO–LUMO gap",value:r===null?"n/a":W(r-e)},{label:"Dipole moment",value:`${o.toFixed(3)} D`},{label:"SCF",value:`${t.converged?"converged":"NOT converged"} in ${t.iterations} iterations, ${Math.round(t.millis)} ms`},{label:"Basis functions",value:String(t.nbf)}];return t.maxForce!==void 0&&n.push({label:"Max force",value:`${t.maxForce.toExponential(2)} Ha/bohr`}),n}const Je=1.6;function Qe(t,e=8,r=6){const n=(t.unrestricted?["alpha","beta"]:["alpha"]).flatMap(u=>{const d=u==="alpha"?t.nalpha:t.nbeta,p=de(t,u),m=p[d-1]?.energy??p[0].energy;return p.filter(h=>h.index>=d-e&&h.index<=d+r-1&&h.energy>m-Je).map(h=>({...h,spin:u}))}),i=n.map(u=>u.energy),l=Math.min(...i),c=Math.max(...i),f=c-l||1;return n.map(u=>({...u,y:.04+.92*((c-u.energy)/f)}))}const U="http://www.w3.org/2000/svg";function et(t,e,r,o){t.replaceChildren();const n=240,i=300;t.setAttribute("viewBox",`0 0 ${n} ${i}`);const l=Qe(e),c=e.unrestricted?2:1,f=70,u=p=>c===1?70:p==="alpha"?40:125,d=new Map;for(const p of l){const m=`${p.spin}:${Math.round(p.energy*k*50)}`;d.set(m,[...d.get(m)??[],p])}for(const p of d.values()){const m=(f-6*(p.length-1))/p.length;p.forEach((h,b)=>{const x=u(h.spin)+b*(m+6),C=h.y*(i-20)+10,v=r?.spin===h.spin&&r.index===h.index,g=document.createElementNS(U,"g");g.classList.add("level"),h.occupied&&g.classList.add("occupied"),v&&g.classList.add("selected");const z=document.createElementNS(U,"rect");z.setAttribute("x",String(x-2)),z.setAttribute("y",String(C-7)),z.setAttribute("width",String(m+4)),z.setAttribute("height","14"),z.setAttribute("class","hit");const w=document.createElementNS(U,"line");w.setAttribute("x1",String(x)),w.setAttribute("x2",String(x+m)),w.setAttribute("y1",String(C)),w.setAttribute("y2",String(C));const G=document.createElementNS(U,"title");if(G.textContent=`${h.label}${c===2?` (${h.spin==="alpha"?"α":"β"})`:""}: ${(h.energy*k).toFixed(2)} eV`,g.append(z,w,G),h.occupied){const B=e.unrestricted?1:2;for(let V=0;V<B;V++){const T=document.createElementNS(U,"text");T.setAttribute("x",String(x+m/2+(B===2?V===0?-5:5:0))),T.setAttribute("y",String(C+4)),T.setAttribute("class","electron"),T.textContent=h.spin==="beta"!=(V===1)?"↓":"↑",g.append(T)}}if(b===p.length-1&&(h.label==="HOMO"||h.label==="LUMO"||v)){const B=document.createElementNS(U,"text");B.setAttribute("x",String(x+m+6)),B.setAttribute("y",String(C+4)),B.setAttribute("class","label"),B.textContent=h.label,g.append(B)}g.addEventListener("click",()=>o({spin:h.spin,index:h.index})),t.append(g)})}if(c===2)for(const[p,m]of[["alpha","α"],["beta","β"]]){const h=document.createElementNS(U,"text");h.setAttribute("x",String(u(p)+f/2)),h.setAttribute("y",String(i-2)),h.setAttribute("class","spin-label"),h.textContent=m,t.append(h)}}const tt=.17;function me(t){return .35+.3*N(t).covalentRadius*D}function X(t){return N(t).color}function nt(t,e){const r=t.map((n,i)=>({center:n.xyz,radius:me(n.symbol),color:X(n.symbol),highlight:i===e})),o=Le(t).map(([n,i])=>({a:t[n].xyz,b:t[i].xyz,radius:tt,colorA:X(t[n].symbol),colorB:X(t[i].symbol)}));return{spheres:r,bonds:o}}function rt(t){const e=[0,1,2].map(o=>t.reduce((n,i)=>n+i.xyz[o],0)/Math.max(1,t.length)),r=t.reduce((o,n)=>Math.max(o,Math.hypot(n.xyz[0]-e[0],n.xyz[1]-e[1],n.xyz[2]-e[2])),0);return{center:e,radius:r}}const y=t=>document.getElementById(t),a={canvas:y("view"),noGpu:y("no-gpu"),preset:y("preset"),xyzBox:y("xyz-box"),xyz:y("xyz"),loadXyz:y("load-xyz"),charge:y("charge"),multiplicity:y("multiplicity"),method:y("method"),basis:y("basis"),relax:y("relax"),resetGeometry:y("reset-geometry"),showCloud:y("show-cloud"),cloudOpacity:y("cloud-opacity"),cloudMin:y("cloud-min"),showOrbital:y("show-orbital"),iso:y("iso"),status:y("status"),readouts:y("readouts"),orbitalName:y("orbital-name"),diagram:document.getElementById("diagram"),saveImage:y("save-image"),playScf:y("play-scf"),scfFrame:y("scf-frame"),scfLabel:y("scf-label")},s={atoms:[],initialAtoms:[],charge:0,multiplicity:1,basis:"6-31g",method:"hf",result:null,orbital:null,frame:null,relaxId:null,hover:null},A=new ke,F=new Te(new Worker(new URL(""+new URL("engine.worker-Dak6_Vu_.js",import.meta.url).href,import.meta.url),{type:"module"}));let S=null;const ot=new URLSearchParams(location.search).has("headless"),ge=Math.max(16,Math.min(160,Number(new URLSearchParams(location.search).get("res"))||112));let O=!0,$=null;function ye(){return{showCloud:a.showCloud.checked,showOrbital:a.showOrbital.checked&&s.orbital!==null,cloudOpacity:40*Number(a.cloudOpacity.value)**2,logMin:Number(a.cloudMin.value),logMax:1,iso:Number(a.iso.value),orbitalOpacity:.55}}function P(t,e="ok"){a.status.textContent=t,a.status.className=`status ${e==="ok"?"":e}`}function Q(t={}){return{atoms:s.atoms.map(e=>({symbol:e.symbol,xyz:e.xyz})),charge:s.charge,multiplicity:s.multiplicity,basis:s.basis,method:s.method,...t}}function H(){if(!S)return;const{spheres:t,bonds:e}=nt(s.atoms,s.hover);S.setAtoms(t,e),O=!0}function ee(){const t=s.result;if(!S||!t)return;const e=s.frame===null?null:t.history[s.frame];try{S.updateVolume(Ze(t,s.orbital,e??t,ge))}catch(r){P(r instanceof Error?r.message:String(r),"error")}O=!0}function be(){S&&O&&(O=!1,it(),S.render(A,ye())),requestAnimationFrame(be)}function it(){const t=Math.min(window.devicePixelRatio||1,2),e=Math.max(1,Math.floor(a.canvas.clientWidth*t)),r=Math.max(1,Math.floor(a.canvas.clientHeight*t));(a.canvas.width!==e||a.canvas.height!==r)&&(a.canvas.width=e,a.canvas.height=r)}function ve(t){return{spin:"alpha",index:t.nalpha-1}}function te(t,e=!1){const r=s.result&&s.result.nmo===t.nmo&&s.result.unrestricted===t.unrestricted;s.result=t,(!r||!s.orbital||s.orbital.index>=t.nmo)&&(s.orbital=ve(t)),e||I(),s.frame=null,a.scfFrame.max=String(Math.max(0,t.history.length-1)),a.scfFrame.value=a.scfFrame.max,a.scfFrame.disabled=t.history.length===0,a.scfLabel.textContent=t.history.length?`iteration ${t.history.length} of ${t.history.length}`:"SCF iterations",a.readouts.replaceChildren(...Ke(t).flatMap(o=>{const n=document.createElement("dt");n.textContent=o.label;const i=document.createElement("dd");return i.textContent=o.value,o.title&&(i.title=o.title),[n,i]})),xe(),ee()}function xe(){const t=s.result;if(t&&(et(a.diagram,t,s.orbital,e=>{s.orbital=e,a.showOrbital.checked=!0,xe(),ee()}),s.orbital)){const e=de(t,s.orbital.spin)[s.orbital.index];a.orbitalName.textContent=`· ${e.label}${t.unrestricted?s.orbital.spin==="alpha"?" α":" β":""}`}}function ne(t){return`${t.energy.toFixed(6)} Ha · ${t.iterations} iterations · ${Math.round(t.millis)} ms`}async function R(t={}){if(s.atoms.length!==0){P("Solving the SCF…","busy");try{const e=await F.computeLatest(Q({history:t.history??!0}));if(!e)return;te(e),P(e.converged?ne(e):`SCF did not converge (${e.iterations} iterations)`,e.converged?"ok":"error")}catch(e){P(e instanceof Error?e.message:String(e),"error")}}}function we(t,e,r){_(),s.atoms=t,s.initialAtoms=t.map(i=>({symbol:i.symbol,xyz:[...i.xyz]})),s.charge=e,s.multiplicity=r,s.result=null,s.orbital=null,a.charge.value=String(e),a.multiplicity.value=String(r),a.xyz.value=Fe(t).split(`
`).slice(2).join(`
`),re(),F.reset();const{center:o,radius:n}=rt(t);A.fit(o,n),H(),R()}function Me(t){const e=Ge(t);for(const r of Array.from(a.basis.options))r.disabled=!e.bases.includes(r.value);e.bases.includes(s.basis)||(s.basis=e.bases[e.bases.length-1],a.basis.value=s.basis),we(ue(e.xyz),e.charge,e.multiplicity)}function re(){const t=s.multiplicity!==1;for(const e of Array.from(a.method.options))e.disabled=t&&e.value!=="hf";t&&s.method!=="hf"&&(s.method="hf",a.method.value="hf")}function _(){s.relaxId!==null&&F.cancelRelax(s.relaxId),s.relaxId=null,a.relax.textContent="Relax geometry"}function at(){if(s.relaxId!==null)return _();a.relax.textContent="Stop",P("Relaxing…","busy"),s.relaxId=F.relax(Q({history:!1}),{onStep:(t,e)=>{s.atoms=s.atoms.map((r,o)=>({symbol:r.symbol,xyz:t.coordinates.slice(3*o,3*o+3)})),H(),te(t),P(`Relaxing: step ${e+1}, E = ${t.energy.toFixed(6)} Ha, max force ${t.maxForce?.toExponential(1)}`,"busy")},onDone:(t,e)=>{s.relaxId=null,a.relax.textContent="Relax geometry",P(t?`Relaxed in ${e} steps · ${ne(s.result)}`:`Stopped after ${e} steps without converging`,t?"ok":"error")},onError:t=>{s.relaxId=null,a.relax.textContent="Relax geometry",P(t,"error")}})}function I(){$!==null&&clearInterval($),$=null,a.playScf.textContent="▶ Play SCF"}function Z(t){const e=s.result;if(!e||e.history.length===0)return;s.frame=t>=e.history.length-1?null:t,a.scfFrame.value=String(t);const r=e.history[t];a.scfLabel.textContent=`iteration ${t+1} of ${e.history.length} · E = ${r.energy.toFixed(6)} · error ${r.error.toExponential(1)}`,ee()}async function st(){if($!==null)return I();_(),P("Recording SCF iterations…","busy");const t=await F.compute(Q({history:!0,warmStart:!1})).catch(r=>(P(r.message,"error"),null));if(!t)return;te(t,!0),P(ne(t));let e=0;a.playScf.textContent="■ Stop",Z(0),$=window.setInterval(()=>{if(e++,e>=t.history.length)return I();Z(e)},450)}function Pe(t){const e=a.canvas.getBoundingClientRect();return[(t.clientX-e.left)/e.width*2-1,1-(t.clientY-e.top)/e.height*2]}function Be(t){const{origin:e,dir:r}=A.ray(...Pe(t));let o=null,n=1/0;return s.atoms.forEach((i,l)=>{const c=Ve(e,r,i.xyz,me(i.symbol)*1.1);c!==null&&c<n&&(n=c,o=l)}),o}let M=null;a.canvas.addEventListener("pointerdown",t=>{a.canvas.setPointerCapture(t.pointerId);const e=Be(t);e!==null?(_(),I(),M={kind:"atom",index:e,planePoint:[...s.atoms[e].xyz],normal:A.basis().forward}):M={kind:"orbit",x:t.clientX,y:t.clientY},a.canvas.classList.add("dragging")});a.canvas.addEventListener("pointermove",t=>{if(!M){const c=Be(t);a.canvas.classList.toggle("over-atom",c!==null),c!==s.hover&&(s.hover=c,H());return}if(M.kind==="orbit"){A.orbit(t.clientX-M.x,t.clientY-M.y),M.x=t.clientX,M.y=t.clientY,O=!0;return}const{origin:e,dir:r}=A.ray(...Pe(t)),o=L(r,M.normal);if(Math.abs(o)<1e-6)return;const n=L(q(M.planePoint,e),M.normal)/o,i=he(e,K(r,n)),l=M.index;s.atoms=s.atoms.map((c,f)=>f===l?{symbol:c.symbol,xyz:i}:c),H(),R({history:!1})});const ze=()=>{M=null,a.canvas.classList.remove("dragging")};a.canvas.addEventListener("pointerup",ze);a.canvas.addEventListener("pointercancel",ze);a.canvas.addEventListener("wheel",t=>{t.preventDefault(),A.zoom(Math.exp(t.deltaY*.001)),O=!0},{passive:!1});window.addEventListener("resize",()=>O=!0);for(const t of fe)a.preset.add(new Option(t.name,t.id));a.preset.add(new Option("Custom (pasted XYZ)","custom"));a.preset.addEventListener("change",()=>{a.preset.value==="custom"?a.xyzBox.open=!0:Me(a.preset.value)});a.loadXyz.addEventListener("click",()=>{try{const t=ue(a.xyz.value);if(t.length===0)throw new Error('No atoms found. Use lines like "O 0 0 0" in Angstrom.');a.preset.value="custom";for(const r of Array.from(a.basis.options))r.disabled=!1;const e=ce({atoms:t,charge:s.charge,multiplicity:1});we(t,s.charge,e%2===0?1:2)}catch(t){P(t instanceof Error?t.message:String(t),"error")}});a.charge.addEventListener("change",()=>{s.charge=Math.round(Number(a.charge.value)),(ce({atoms:s.atoms,charge:s.charge})+s.multiplicity-1)%2!==0&&(s.multiplicity=s.multiplicity===1?2:s.multiplicity-1),a.multiplicity.value=String(s.multiplicity),re(),R()});a.multiplicity.addEventListener("change",()=>{s.multiplicity=Math.max(1,Math.round(Number(a.multiplicity.value))),re(),R()});a.method.addEventListener("change",()=>{s.method=a.method.value,R()});a.basis.addEventListener("change",()=>{s.basis=a.basis.value,R()});a.relax.addEventListener("click",at);a.resetGeometry.addEventListener("click",()=>{_(),s.atoms=s.initialAtoms.map(t=>({symbol:t.symbol,xyz:[...t.xyz]})),H(),R()});for(const t of[a.showCloud,a.cloudOpacity,a.cloudMin,a.iso])t.addEventListener("input",()=>O=!0);a.showOrbital.addEventListener("input",()=>{a.showOrbital.checked&&!s.orbital&&s.result&&(s.orbital=ve(s.result)),O=!0});a.playScf.addEventListener("click",st);a.saveImage.addEventListener("click",async()=>{const t=Math.min(window.devicePixelRatio||1,2),e=document.createElement("a");e.href=await Ee(Math.round(a.canvas.clientWidth*t),Math.round(a.canvas.clientHeight*t)),e.download=`${a.preset.value}-${s.method}-${s.basis}.png`,e.click()});a.scfFrame.addEventListener("input",()=>{I(),Z(Number(a.scfFrame.value))});async function lt(){a.method.value=s.method,a.basis.value=s.basis;try{S=await J.create(ot?null:a.canvas),S.device.lost.then(t=>P(`The GPU device was lost: ${t.message}`,"error"))}catch(t){a.noGpu.hidden=!1,a.noGpu.textContent=`${t instanceof Error?t.message:String(t)} The calculations still run and the numbers on the right are real, but the 3D view needs a browser with WebGPU (recent Chrome, Edge, or Safari).`}a.preset.value="h2o",await F.ready,Me("h2o"),requestAnimationFrame(be)}async function Ee(t=a.canvas.clientWidth||800,e=a.canvas.clientHeight||600){if(!S)throw new Error("WebGPU is not available.");A.aspect=t/e;const r=await S.snapshot(A,ye(),t,e),o=document.createElement("canvas");return o.width=t,o.height=e,o.getContext("2d").putImageData(new ImageData(r,t,e),0,0),o.toDataURL("image/png")}Object.assign(window,{__qm:{state:s,engine:F,camera:A,snapshot:Ee,volumeResolution:ge,get renderer(){return S}}});lt();
