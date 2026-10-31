(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const r of document.querySelectorAll('link[rel="modulepreload"]'))o(r);new MutationObserver(r=>{for(const i of r)if(i.type==="childList")for(const l of i.addedNodes)l.tagName==="LINK"&&l.rel==="modulepreload"&&o(l)}).observe(document,{childList:!0,subtree:!0});function n(r){const i={};return r.integrity&&(i.integrity=r.integrity),r.referrerPolicy&&(i.referrerPolicy=r.referrerPolicy),r.crossOrigin==="use-credentials"?i.credentials="include":r.crossOrigin==="anonymous"?i.credentials="omit":i.credentials="same-origin",i}function o(r){if(r.ep)return;r.ep=!0;const i=n(r);fetch(r.href,i)}})();const j=1/.52917721092,q=27.211386245988,$e=2.541746473,Ie=[["H",1.008,.31,[.95,.95,.95]],["He",4.0026,.28,[.85,1,1]],["Li",6.94,1.28,[.8,.5,1]],["Be",9.0122,.96,[.76,1,0]],["B",10.81,.84,[1,.71,.71]],["C",12.011,.76,[.38,.38,.38]],["N",14.007,.71,[.19,.31,.97]],["O",15.999,.66,[1,.13,.13]],["F",18.998,.57,[.56,.88,.31]],["Ne",20.18,.58,[.7,.89,.96]]],He=Ie.map(([t,e,n,o],r)=>({symbol:t,z:r+1,mass:e,covalentRadius:n,color:o})),_e=new Map(He.map(t=>[t.symbol.toLowerCase(),t]));function C(t){const e=_e.get(t.toLowerCase());if(!e)throw new Error(`Unsupported element: ${t} (only H to Ne are available)`);return e}function Ve(t){return C(t.symbol).z}function we(t){return t.atoms.reduce((e,n)=>e+Ve(n),0)-t.charge}function De(t,e=1.25){const n=[];for(let o=0;o<t.length;o++)for(let r=0;r<o;r++){const[i,l,c]=t[o].xyz,[h,u,d]=t[r].xyz;Math.hypot(i-h,l-u,c-d)/j<e*(C(t[o].symbol).covalentRadius+C(t[r].symbol).covalentRadius)&&n.push([r,o])}return n}function Me(t){const e=[];for(const n of t.split(/\r?\n/)){const o=n.trim().split(/\s+/);if(o.length<4)continue;const[r,...i]=o,l=i.slice(0,3).map(Number);if(!/^[A-Za-z]{1,2}$/.test(r)||l.some(h=>!Number.isFinite(h)))continue;const c=r[0].toUpperCase()+r.slice(1).toLowerCase();C(c),e.push({symbol:c,xyz:l.map(h=>h*j)})}return e}function qe(t,e=""){const n=t.map(o=>`${o.symbol.padEnd(2)} ${o.xyz.map(r=>(r/j).toFixed(6).padStart(12)).join(" ")}`);return`${t.length}
${e}
${n.join(`
`)}
`}const O=["sto-3g","6-31g","6-31g*"];function je(){const t=[];for(let e=0;e<6;e++){const n=e*Math.PI/3;t.push(`C ${(1.397*Math.cos(n)).toFixed(4)} ${(1.397*Math.sin(n)).toFixed(4)} 0`)}for(let e=0;e<6;e++){const n=e*Math.PI/3;t.push(`H ${(2.481*Math.cos(n)).toFixed(4)} ${(2.481*Math.sin(n)).toFixed(4)} 0`)}return t.join(`
`)}const re=[{id:"h2",name:"Hydrogen H₂",xyz:`H 0 0 -0.3705
H 0 0 0.3705`,charge:0,multiplicity:1,bases:O},{id:"h2o",name:"Water H₂O",xyz:`O 0 0 0.1173
H 0 0.7572 -0.4692
H 0 -0.7572 -0.4692`,charge:0,multiplicity:1,bases:O},{id:"nh3",name:"Ammonia NH₃",xyz:`N 0 0 0.1162
H 0 0.9377 -0.2711
H 0.8121 -0.4689 -0.2711
H -0.8121 -0.4689 -0.2711`,charge:0,multiplicity:1,bases:O},{id:"ch4",name:"Methane CH₄",xyz:`C 0 0 0
H 0.6291 0.6291 0.6291
H -0.6291 -0.6291 0.6291
H -0.6291 0.6291 -0.6291
H 0.6291 -0.6291 -0.6291`,charge:0,multiplicity:1,bases:O},{id:"hf",name:"Hydrogen fluoride HF",xyz:`F 0 0 0.0917
H 0 0 -0.8253`,charge:0,multiplicity:1,bases:O},{id:"n2",name:"Nitrogen N₂",xyz:`N 0 0 -0.5488
N 0 0 0.5488`,charge:0,multiplicity:1,bases:O},{id:"co2",name:"Carbon dioxide CO₂",xyz:`C 0 0 0
O 0 0 1.1621
O 0 0 -1.1621`,charge:0,multiplicity:1,bases:O},{id:"ch2o",name:"Formaldehyde CH₂O",xyz:`C 0 0 -0.5285
O 0 0 0.6745
H 0 0.9363 -1.1079
H 0 -0.9363 -1.1079`,charge:0,multiplicity:1,bases:O},{id:"c2h4",name:"Ethylene C₂H₄",xyz:`C 0 0 0.6695
C 0 0 -0.6695
H 0 0.9288 1.2321
H 0 -0.9288 1.2321
H 0 0.9288 -1.2321
H 0 -0.9288 -1.2321`,charge:0,multiplicity:1,bases:O},{id:"c6h6",name:"Benzene C₆H₆",xyz:je(),charge:0,multiplicity:1,bases:["sto-3g"]},{id:"oh",name:"Hydroxyl radical OH",xyz:`O 0 0 0
H 0 0 0.9697`,charge:0,multiplicity:2,bases:O},{id:"o2",name:"Oxygen O₂ (triplet)",xyz:`O 0 0 -0.6038
O 0 0 0.6038`,charge:0,multiplicity:3,bases:O}];function We(t){const e=re.find(n=>n.id===t);if(!e)throw new Error(`Unknown preset ${t}`);return e}const Ye={h2:"The simplest molecule: two hydrogen atoms sharing their two electrons in one bond.",h2o:"Oxygen pulls the shared electrons toward itself and leaves the hydrogens slightly positive. That imbalance is what makes water such a good solvent.",nh3:"Nitrogen bonds three hydrogens and keeps one pair of electrons to itself, which pushes the molecule into a pyramid.",ch4:"Carbon shares its four outer electrons with four hydrogens, which point to the corners of a tetrahedron.",hf:"Fluorine pulls on electrons harder than any other element, so this bond is very lopsided.",n2:"A triple bond, one of the strongest in chemistry. It is why the nitrogen that makes up most of the air hardly reacts.",co2:"Two double bonds in a straight line. Each bond is lopsided, but they point in opposite directions and cancel out.",ch2o:"A carbon–oxygen double bond, the group at the heart of sugars and many other molecules of life.",c2h4:"A carbon–carbon double bond. Show the HOMO orbital to see the extra bond lying above and below the molecule.",c6h6:"Six electrons shared around the whole ring instead of sitting in fixed bonds. It uses the quickest settings to stay smooth.",oh:"A radical: one electron has no partner, which makes it extremely reactive.",o2:"The oxygen we breathe has two unpaired electrons, which is why liquid oxygen sticks to a magnet."};function Ee(t,e="alpha"){const n=e==="beta"&&t.unrestricted?t.epsB:t.epsA,o=e==="beta"?t.nbeta:t.nalpha;return n.map((r,i)=>{let l="";return i===o-1?l="HOMO":i===o?l="LUMO":i<o?l=`HOMO−${o-1-i}`:l=`LUMO+${i-o}`,{index:i,energy:r,occupied:i<o,label:l}})}function Xe(t){let e=t.epsA[t.nalpha-1]??-1/0,n=t.epsA[t.nalpha]??1/0;return t.unrestricted&&(t.nbeta>0&&(e=Math.max(e,t.epsB[t.nbeta-1])),n=Math.min(n,t.epsB[t.nbeta]??1/0)),{homo:e,lumo:Number.isFinite(n)?n:null}}class Ke{worker;nextId=1;pending=new Map;relaxes=new Map;busy=!1;queued=null;ready;constructor(e){this.worker=e;let n;this.ready=new Promise(o=>n=o),e.onmessage=o=>{const r=o.data;if(r.type==="ready")return n();if(r.type==="relax-step"||r.type==="relax-done"){const l=this.relaxes.get(r.id);if(!l)return;r.type==="relax-step"?l.onStep(r.output,r.step):(this.relaxes.delete(r.id),l.onDone(r.converged,r.steps));return}if(r.type==="error"&&this.relaxes.has(r.id)){const l=this.relaxes.get(r.id);this.relaxes.delete(r.id),l.onError(r.message);return}const i=this.pending.get(r.id);i&&(this.pending.delete(r.id),r.type==="result"?i.resolve(r.output):r.type==="error"&&i.reject(new Error(r.message)))}}send(e){this.worker.postMessage(e)}compute(e){const n=this.nextId++;return new Promise((o,r)=>{this.pending.set(n,{resolve:o,reject:r}),this.send({id:n,type:"compute",input:e})})}computeLatest(e){return this.busy?(this.queued?.resolve(null),new Promise(n=>this.queued={input:e,resolve:n})):(this.busy=!0,this.compute(e).finally(()=>{this.busy=!1;const n=this.queued;this.queued=null,n&&this.computeLatest(n.input).then(n.resolve,()=>n.resolve(null))}))}relax(e,n,o=60){const r=this.nextId++;return this.relaxes.set(r,n),this.send({id:r,type:"relax",input:e,maxSteps:o}),r}cancelRelax(e){this.relaxes.delete(e),this.send({id:e,type:"cancel"})}reset(){this.send({id:this.nextId++,type:"reset"})}}const W=(t,e)=>[t[0]-e[0],t[1]-e[1],t[2]-e[2]],Be=(t,e)=>[t[0]+e[0],t[1]+e[1],t[2]+e[2]],oe=(t,e)=>[t[0]*e,t[1]*e,t[2]*e],T=(t,e)=>t[0]*e[0]+t[1]*e[1]+t[2]*e[2],pe=(t,e)=>[t[1]*e[2]-t[2]*e[1],t[2]*e[0]-t[0]*e[2],t[0]*e[1]-t[1]*e[0]],Ze=t=>Math.hypot(t[0],t[1],t[2]),te=t=>oe(t,1/(Ze(t)||1));function Je(t,e){const n=new Float32Array(16);for(let o=0;o<4;o++)for(let r=0;r<4;r++){let i=0;for(let l=0;l<4;l++)i+=t[l*4+r]*e[o*4+l];n[o*4+r]=i}return n}function Qe(t,e,n,o){const r=1/Math.tan(t/2),i=new Float32Array(16);return i[0]=r/e,i[5]=r,i[10]=o/(n-o),i[11]=-1,i[14]=o*n/(n-o),i}function et(t,e,n){const o=te(W(t,e)),r=te(pe(n,o)),i=pe(o,r);return new Float32Array([r[0],i[0],o[0],0,r[1],i[1],o[1],0,r[2],i[2],o[2],0,-T(r,t),-T(i,t),-T(o,t),1])}function Pe(t){const e=new Float64Array(16),n=Array.from(t),o=Array.from({length:4},(r,i)=>[0,1,2,3].map(l=>n[l*4+i]).concat([0,1,2,3].map(l=>l===i?1:0)));for(let r=0;r<4;r++){let i=r;for(let c=r+1;c<4;c++)Math.abs(o[c][r])>Math.abs(o[i][r])&&(i=c);[o[r],o[i]]=[o[i],o[r]];const l=o[r][r];if(Math.abs(l)<1e-30)throw new Error("singular matrix");for(let c=0;c<8;c++)o[r][c]/=l;for(let c=0;c<4;c++){if(c===r)continue;const h=o[c][r];for(let u=0;u<8;u++)o[c][u]-=h*o[r][u]}}for(let r=0;r<4;r++)for(let i=0;i<4;i++)e[i*4+r]=o[r][4+i];return new Float32Array(e)}function me(t,e){const n=t[0]*e[0]+t[4]*e[1]+t[8]*e[2]+t[12],o=t[1]*e[0]+t[5]*e[1]+t[9]*e[2]+t[13],r=t[2]*e[0]+t[6]*e[1]+t[10]*e[2]+t[14],i=t[3]*e[0]+t[7]*e[1]+t[11]*e[2]+t[15];return[n/i,o/i,r/i]}function tt(t,e,n,o){const r=W(t,n),i=T(r,e),l=T(r,r)-o*o,c=i*i-l;if(c<0)return null;const h=-i-Math.sqrt(c);return h>0?h:null}class nt{target=[0,0,0];distance=15;yaw=.6;pitch=.35;fovy=40*Math.PI/180;aspect=1;get eye(){const e=Math.cos(this.pitch);return Be(this.target,oe([e*Math.sin(this.yaw),Math.sin(this.pitch),e*Math.cos(this.yaw)],this.distance))}get view(){return et(this.eye,this.target,[0,1,0])}get projection(){return Qe(this.fovy,this.aspect,Math.max(.05,this.distance*.02),this.distance*10+50)}get viewProjection(){return Je(this.projection,this.view)}orbit(e,n){this.yaw-=e*.008,this.pitch=Math.max(-1.5,Math.min(1.5,this.pitch+n*.008))}zoom(e){this.distance=Math.max(2,Math.min(200,this.distance*e))}fit(e,n){this.target=e,this.distance=Math.max(6,(n+2.5)/Math.sin(this.fovy/2))}ray(e,n){const o=Pe(this.viewProjection),r=me(o,[e,n,0]),i=me(o,[e,n,1]);return{origin:r,dir:te(W(i,r))}}basis(){const e=this.view;return{right:[e[0],e[4],e[8]],up:[e[1],e[5],e[9]],forward:[-e[2],-e[6],-e[10]]}}}const rt=`// Ball-and-stick atoms: spheres are ray-cast impostors, bonds are instanced cylinders.

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
`,ot=`// Ray-marches the density cloud and the orbital isosurfaces, stopping at the atoms' depth.

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
`,it=`// Evaluates the electron density and one molecular orbital on a 3D grid.
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
`,ge="depth24plus",J=18;function st(){const t=[];for(let e=0;e<J;e++){const n=e/J*2*Math.PI,o=(e+1)/J*2*Math.PI,r=(i,l)=>t.push(Math.cos(i),Math.sin(i),l);r(n,0),r(o,0),r(o,1),r(n,0),r(o,1),r(n,1)}return new Float32Array(t)}class ie{device;context;format;depth=null;cameraBuffer;rayBuffer;paramsBuffer;normsBuffer;ringBuffer;ringCount;sphereBuffer=null;sphereCount=0;bondBuffer=null;bondCount=0;spherePipeline;bondPipeline;rayPipeline;volumePipeline;atomsBindGroup;sampler;volume=null;box=null;frame=0;adapter=null;background={r:.035,g:.045,b:.07,a:1};static async create(e){if(!navigator.gpu)throw new Error("WebGPU is not available in this browser.");const n=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!n)throw new Error("No WebGPU adapter was found.");const o=await n.requestDevice(),r=e?e.getContext("webgpu"):null;if(e&&!r)throw new Error("Could not create a WebGPU canvas context.");const i=new ie(o,r,navigator.gpu.getPreferredCanvasFormat());return i.adapter=n,i}constructor(e,n,o){this.device=e,this.context=n,this.format=o,n?.configure({device:e,format:o,alphaMode:"opaque"});const r=GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST;this.cameraBuffer=e.createBuffer({size:128,usage:r}),this.rayBuffer=e.createBuffer({size:144,usage:r}),this.paramsBuffer=e.createBuffer({size:48,usage:r}),this.normsBuffer=e.createBuffer({size:40,usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});const i=st();this.ringCount=i.length/3,this.ringBuffer=e.createBuffer({size:i.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),e.queue.writeBuffer(this.ringBuffer,0,i),this.sampler=e.createSampler({magFilter:"linear",minFilter:"linear"});const l=e.createShaderModule({code:rt,label:"atoms"}),c=e.createPipelineLayout({bindGroupLayouts:[this.cameraLayout()]}),h={format:ge,depthWriteEnabled:!0,depthCompare:"less"},u=(g,f)=>({shaderLocation:g,offset:f,format:"float32x3"}),d=(g,f)=>({shaderLocation:g,offset:f,format:"float32"});this.spherePipeline=e.createRenderPipeline({label:"spheres",layout:c,vertex:{module:l,entryPoint:"sphereVs",buffers:[{arrayStride:32,stepMode:"instance",attributes:[u(0,0),d(1,12),u(2,16),d(3,28)]}]},fragment:{module:l,entryPoint:"sphereFs",targets:[{format:o}]},primitive:{topology:"triangle-strip"},depthStencil:h}),this.bondPipeline=e.createRenderPipeline({label:"bonds",layout:c,vertex:{module:l,entryPoint:"bondVs",buffers:[{arrayStride:52,stepMode:"instance",attributes:[u(0,0),d(1,12),u(2,16),u(3,28),u(4,40)]},{arrayStride:12,stepMode:"vertex",attributes:[u(5,0)]}]},fragment:{module:l,entryPoint:"bondFs",targets:[{format:o}]},primitive:{topology:"triangle-list",cullMode:"back"},depthStencil:h}),this.atomsBindGroup=e.createBindGroup({layout:this.cameraLayout(),entries:[{binding:0,resource:{buffer:this.cameraBuffer}}]});const m=e.createShaderModule({code:ot,label:"raymarch"});this.rayPipeline=e.createRenderPipeline({label:"raymarch",layout:"auto",vertex:{module:m,entryPoint:"vs"},fragment:{module:m,entryPoint:"fs",targets:[{format:o,blend:{color:{srcFactor:"one",dstFactor:"one-minus-src-alpha"},alpha:{srcFactor:"one",dstFactor:"one-minus-src-alpha"}}}]}}),this.volumePipeline=e.createComputePipeline({label:"volume",layout:"auto",compute:{module:e.createShaderModule({code:it,label:"volume"}),entryPoint:"main"}})}cameraLayoutCache=null;cameraLayout(){return this.cameraLayoutCache??=this.device.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.VERTEX|GPUShaderStage.FRAGMENT,buffer:{type:"uniform"}}]}),this.cameraLayoutCache}resize(e,n){this.depth&&this.depth.width===e&&this.depth.height===n||(this.depth?.destroy(),this.depth=this.device.createTexture({size:[e,n],format:ge,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.TEXTURE_BINDING}))}setAtoms(e,n){const o=new Float32Array(Math.max(1,e.length)*8);e.forEach((i,l)=>o.set([...i.center,i.radius,...i.color,i.highlight?1:0],8*l)),(!this.sphereBuffer||this.sphereBuffer.size<o.byteLength)&&(this.sphereBuffer?.destroy(),this.sphereBuffer=this.device.createBuffer({size:o.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST})),this.device.queue.writeBuffer(this.sphereBuffer,0,o),this.sphereCount=e.length;const r=new Float32Array(Math.max(1,n.length)*13);n.forEach((i,l)=>r.set([...i.a,i.radius,...i.b,...i.colorA,...i.colorB],13*l)),(!this.bondBuffer||this.bondBuffer.size<r.byteLength)&&(this.bondBuffer?.destroy(),this.bondBuffer=this.device.createBuffer({size:r.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST})),this.device.queue.writeBuffer(this.bondBuffer,0,r),this.bondCount=n.length}updateVolume(e){const{device:n}=this,[o,r,i]=e.box.dims;(!this.volume||this.volume.width!==o||this.volume.height!==r||this.volume.depthOrArrayLayers!==i)&&(this.volume?.destroy(),this.volume=n.createTexture({size:[o,r,i],dimension:"3d",format:"rgba16float",usage:GPUTextureUsage.STORAGE_BINDING|GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_SRC})),this.box=e.box;const l=new ArrayBuffer(48),c=new Float32Array(l),h=new Uint32Array(l);c.set(e.box.min,0),h[3]=e.nbf,c.set(e.box.max,4),h[7]=e.nshell,h.set(e.box.dims,8),h[11]=e.nocc,n.queue.writeBuffer(this.paramsBuffer,0,l),n.queue.writeBuffer(this.normsBuffer,0,e.norms);const u=y=>{const x=n.createBuffer({size:Math.max(16,y.byteLength),usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST});return n.queue.writeBuffer(x,0,y),x},d=[u(e.shells),u(e.prims),u(e.coefs),u(e.weights.length?e.weights:new Float32Array(4))],m=n.createBindGroup({layout:this.volumePipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.paramsBuffer}},{binding:1,resource:{buffer:d[0]}},{binding:2,resource:{buffer:d[1]}},{binding:3,resource:{buffer:this.normsBuffer}},{binding:4,resource:{buffer:d[2]}},{binding:5,resource:{buffer:d[3]}},{binding:6,resource:this.volume.createView()}]}),g=n.createCommandEncoder(),f=g.beginComputePass();f.setPipeline(this.volumePipeline),f.setBindGroup(0,m),f.dispatchWorkgroups(Math.ceil(o/4),Math.ceil(r/4),Math.ceil(i/4)),f.end(),n.queue.submit([g.finish()]);for(const y of d)y.destroy()}render(e,n,o){const{device:r}=this,i=o??this.context?.getCurrentTexture();if(!i)return;this.resize(i.width,i.height),e.aspect=i.width/i.height;const l=new Float32Array(32);l.set(e.view,0),l.set(e.projection,16),r.queue.writeBuffer(this.cameraBuffer,0,l);const c=r.createCommandEncoder(),h=c.beginRenderPass({colorAttachments:[{view:i.createView(),clearValue:this.background,loadOp:"clear",storeOp:"store"}],depthStencilAttachment:{view:this.depth.createView(),depthClearValue:1,depthLoadOp:"clear",depthStoreOp:"store"}});if(this.bondCount>0&&this.bondBuffer&&(h.setPipeline(this.bondPipeline),h.setBindGroup(0,this.atomsBindGroup),h.setVertexBuffer(0,this.bondBuffer),h.setVertexBuffer(1,this.ringBuffer),h.draw(this.ringCount,this.bondCount)),this.sphereCount>0&&this.sphereBuffer&&(h.setPipeline(this.spherePipeline),h.setBindGroup(0,this.atomsBindGroup),h.setVertexBuffer(0,this.sphereBuffer),h.draw(4,this.sphereCount)),h.end(),this.volume&&this.box&&(n.showCloud||n.showOrbital)){const u=this.box,d=new Float32Array(36);d.set(Pe(e.viewProjection),0),d.set(e.eye,16),d[19]=n.cloudOpacity,d.set(u.min,20),d[23]=n.logMin,d.set(u.max,24),d[27]=n.logMax,d.set(u.max.map((f,y)=>(f-u.min[y])/Math.max(1,u.dims[y]-1)),28),d[31]=n.iso,d[32]=n.showCloud?1:0,d[33]=n.showOrbital?1:0,d[34]=n.orbitalOpacity,d[35]=this.frame++%64,r.queue.writeBuffer(this.rayBuffer,0,d);const m=r.createBindGroup({layout:this.rayPipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.rayBuffer}},{binding:1,resource:this.volume.createView()},{binding:2,resource:this.sampler},{binding:3,resource:this.depth.createView()}]}),g=c.beginRenderPass({colorAttachments:[{view:i.createView(),loadOp:"load",storeOp:"store"}]});g.setPipeline(this.rayPipeline),g.setBindGroup(0,m),g.draw(3),g.end()}r.queue.submit([c.finish()])}async snapshot(e,n,o,r){const i=this.device.createTexture({size:[o,r],format:this.format,usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});this.render(e,n,i);const l=Math.ceil(o*4/256)*256,c=this.device.createBuffer({size:l*r,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),h=this.device.createCommandEncoder();h.copyTextureToBuffer({texture:i},{buffer:c,bytesPerRow:l},[o,r]),this.device.queue.submit([h.finish()]),await c.mapAsync(GPUMapMode.READ);const u=new Uint8Array(c.getMappedRange()),d=new Uint8ClampedArray(o*r*4),m=this.format.startsWith("bgra");for(let g=0;g<r;g++)for(let f=0;f<o;f++){const y=g*l+4*f,x=4*(g*o+f);d[x]=u[y+(m?2:0)],d[x+1]=u[y+1],d[x+2]=u[y+(m?0:2)],d[x+3]=255}return c.unmap(),c.destroy(),i.destroy(),d}async readVolume(){if(!this.volume||!this.box)return null;const[e,n,o]=this.box.dims,r=Math.ceil(e*8/256)*256,i=this.device.createBuffer({size:r*n*o,usage:GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ}),l=this.device.createCommandEncoder();l.copyTextureToBuffer({texture:this.volume},{buffer:i,bytesPerRow:r,rowsPerImage:n},[e,n,o]),this.device.queue.submit([l.finish()]),await i.mapAsync(GPUMapMode.READ);const c=new Uint16Array(i.getMappedRange()),h=new Float32Array(e*n*o*2);for(let u=0;u<o;u++)for(let d=0;d<n;d++)for(let m=0;m<e;m++){const g=(u*n+d)*r/2+m*4,f=2*((u*n+d)*e+m);h[f]=be(c[g]),h[f+1]=be(c[g+1])}return i.unmap(),i.destroy(),{dims:[e,n,o],data:h}}}function be(t){const e=t&32768?-1:1,n=t>>10&31,o=t&1023;return n===0?e*o*2**-24:n===31?o?NaN:e*(1/0):e*(1+o/1024)*2**(n-15)}const ye=128,at=5;function lt(t,e=128,n=at){const o=[1/0,1/0,1/0],r=[-1/0,-1/0,-1/0];for(let h=0;h<t.length;h+=3)for(let u=0;u<3;u++)o[u]=Math.min(o[u],t[h+u]-n),r[u]=Math.max(r[u],t[h+u]+n);const i=r.map((h,u)=>h-o[u]),l=Math.max(...i)/(e-1),c=i.map(h=>Math.min(e,Math.ceil(h/l-1e-9)+1));for(let h=0;h<3;h++){const u=(o[h]+r[h])/2,d=(c[h]-1)*l/2;o[h]=u-d,r[h]=u+d}return{min:o,max:r,dims:c}}function ct(t,e,n=t,o=128){const{nbf:r,nmo:i}=t;if(r>ye)throw new Error(`The renderer supports at most ${ye} basis functions (this system has ${r}).`);const{shells:l,prims:c,componentNorms:h}=t.basis,u=l.length/8,d=new ArrayBuffer(32*u),m=new Float32Array(d),g=new Uint32Array(d);for(let v=0;v<u;v++){const b=8*v,P=l[b+4],M=l[b+5];let N=1/0;for(let B=0;B<P;B++)N=Math.min(N,c[2*(M+B)]);m[b]=l[b],m[b+1]=l[b+1],m[b+2]=l[b+2],g[b+3]=l[b+3],g[b+4]=P,g[b+5]=M,g[b+6]=l[b+6],m[b+7]=N}const f=t.unrestricted?n.cb:n.ca,y=[];for(let v=0;v<t.nalpha;v++)y.push({c:n.ca,k:v,w:t.unrestricted?1:2});if(t.unrestricted)for(let v=0;v<t.nbeta;v++)y.push({c:f,k:v,w:1});const x=y.length,L=new Float32Array((x+1)*r);if(y.forEach(({c:v,k:b},P)=>{for(let M=0;M<r;M++)L[P*r+M]=v[M*i+b]}),e){const v=e.spin==="beta"?f:n.ca;for(let b=0;b<r;b++)L[x*r+b]=v[b*i+e.index]}return{nbf:r,nshell:u,nocc:x,shells:d,prims:Float32Array.from(c),norms:Float32Array.from(h.flat()),coefs:L,weights:Float32Array.from(y.map(v=>v.w)),box:lt(t.coordinates,o)}}const Q=t=>`${(t*q).toFixed(2)} eV`;function ut(t){const{homo:e,lumo:n}=Xe(t),o=Math.hypot(...t.dipole)*$e,r=[{label:"Total energy",value:`${t.energy.toFixed(6)} Ha`,title:`${(t.energy*q).toFixed(3)} eV. The energy of the electrons plus the repulsion between the nuclei; lower is more stable (1 Ha = 27.211 eV).`},{label:"HOMO",value:Q(e),title:"Energy of the highest occupied orbital. Its negative is roughly the energy needed to remove an electron."},{label:"LUMO",value:n===null?"none":Q(n),title:"Energy of the lowest empty orbital, where an added electron would go."},{label:"HOMO–LUMO gap",value:n===null?"n/a":Q(n-e),title:"Roughly how much energy it takes to excite the molecule."},{label:"Dipole moment",value:`${o.toFixed(3)} D`,title:"How unevenly the charge is spread, in debye. Water is about 1.85 D in experiment."},{label:"SCF",value:`${t.converged?"converged":"NOT converged"} in ${t.iterations} iterations, ${Math.round(t.millis)} ms`,title:"The self-consistent field loop: guess the orbitals, compute the field they make, solve again, until nothing changes."},{label:"Basis functions",value:String(t.nbf),title:"How many atom-centred functions the orbitals are built from. Bigger basis sets are more accurate and slower."}];return t.maxForce!==void 0&&r.push({label:"Max force",value:`${t.maxForce.toExponential(2)} Ha/bohr`}),r}const ht={H:"hydrogen",He:"helium",Li:"lithium",Be:"beryllium",B:"boron",C:"carbon",N:"nitrogen",O:"oxygen",F:"fluorine",Ne:"neon"},dt=([t,e,n])=>`rgb(${Math.round(t*255)}, ${Math.round(e*255)}, ${Math.round(n*255)})`;function ft(t){const e=[...new Set(t.map(n=>C(n.symbol).symbol))];return e.sort((n,o)=>C(o).z-C(n).z),e.map(n=>({symbol:n,name:ht[n]??n,color:dt(C(n).color)}))}function pt(t){return t==="HOMO"?"the highest-energy orbital that holds electrons":t==="LUMO"?"the lowest-energy empty orbital":t.startsWith("HOMO")?"an occupied orbital below the HOMO":"an empty orbital above the LUMO"}const mt=1.6;function gt(t,e=8,n=6){const r=(t.unrestricted?["alpha","beta"]:["alpha"]).flatMap(u=>{const d=u==="alpha"?t.nalpha:t.nbeta,m=Ee(t,u),g=m[d-1]?.energy??m[0].energy;return m.filter(f=>f.index>=d-e&&f.index<=d+n-1&&f.energy>g-mt).map(f=>({...f,spin:u}))}),i=r.map(u=>u.energy),l=Math.min(...i),c=Math.max(...i),h=c-l||1;return r.map(u=>({...u,y:.04+.92*((c-u.energy)/h)}))}const U="http://www.w3.org/2000/svg";function bt(t,e,n,o){t.replaceChildren();const r=240,i=300;t.setAttribute("viewBox",`0 0 ${r} ${i}`);const l=gt(e),c=e.unrestricted?2:1,h=70,u=m=>c===1?70:m==="alpha"?40:125,d=new Map;for(const m of l){const g=`${m.spin}:${Math.round(m.energy*q*50)}`;d.set(g,[...d.get(g)??[],m])}for(const m of d.values()){const g=(h-6*(m.length-1))/m.length;m.forEach((f,y)=>{const x=u(f.spin)+y*(g+6),L=f.y*(i-20)+10,v=n?.spin===f.spin&&n.index===f.index,b=document.createElementNS(U,"g");b.classList.add("level"),f.occupied&&b.classList.add("occupied"),v&&b.classList.add("selected");const P=document.createElementNS(U,"rect");P.setAttribute("x",String(x-2)),P.setAttribute("y",String(L-7)),P.setAttribute("width",String(g+4)),P.setAttribute("height","14"),P.setAttribute("class","hit");const M=document.createElementNS(U,"line");M.setAttribute("x1",String(x)),M.setAttribute("x2",String(x+g)),M.setAttribute("y1",String(L)),M.setAttribute("y2",String(L));const N=document.createElementNS(U,"title");if(N.textContent=`${f.label}${c===2?` (${f.spin==="alpha"?"α":"β"})`:""}: ${(f.energy*q).toFixed(2)} eV`,b.append(P,M,N),f.occupied){const B=e.unrestricted?1:2;for(let D=0;D<B;D++){const G=document.createElementNS(U,"text");G.setAttribute("x",String(x+g/2+(B===2?D===0?-5:5:0))),G.setAttribute("y",String(L+4)),G.setAttribute("class","electron"),G.textContent=f.spin==="beta"!=(D===1)?"↓":"↑",b.append(G)}}if(y===m.length-1&&(f.label==="HOMO"||f.label==="LUMO"||v)){const B=document.createElementNS(U,"text");B.setAttribute("x",String(x+g+6)),B.setAttribute("y",String(L+4)),B.setAttribute("class","label"),B.textContent=f.label,b.append(B)}b.addEventListener("click",()=>o({spin:f.spin,index:f.index})),t.append(b)})}if(c===2)for(const[m,g]of[["alpha","α"],["beta","β"]]){const f=document.createElementNS(U,"text");f.setAttribute("x",String(u(m)+h/2)),f.setAttribute("y",String(i-2)),f.setAttribute("class","spin-label"),f.textContent=g,t.append(f)}}const yt=.17;function Oe(t){return .35+.3*C(t).covalentRadius*j}function ee(t){return C(t).color}function vt(t,e){const n=t.map((r,i)=>({center:r.xyz,radius:Oe(r.symbol),color:ee(r.symbol),highlight:i===e})),o=De(t).map(([r,i])=>({a:t[r].xyz,b:t[i].xyz,radius:yt,colorA:ee(t[r].symbol),colorB:ee(t[i].symbol)}));return{spheres:n,bonds:o}}function xt(t){const e=[0,1,2].map(o=>t.reduce((r,i)=>r+i.xyz[o],0)/Math.max(1,t.length)),n=t.reduce((o,r)=>Math.max(o,Math.hypot(r.xyz[0]-e[0],r.xyz[1]-e[1],r.xyz[2]-e[2])),0);return{center:e,radius:n}}const wt=[{until:"drag",html:"<b>Try it:</b> drag one of the atoms away from the others. Every frame is a fresh quantum calculation."},{until:"relax",html:"The energy on the right went up: the molecule is strained. Press <b>Relax geometry</b> and watch the atoms settle back."},{until:"orbital",html:"Now look inside. Click a line under <b>Energy levels</b> on the right to see where that pair of electrons lives."},{until:"preset",html:"That shape is an orbital: blue and orange are the two signs of its wave. Next, pick another molecule from the list on the left."}];function Mt(t){return wt.find(e=>!t.has(e.until))??null}const p=t=>document.getElementById(t),s={canvas:p("view"),noGpu:p("no-gpu"),preset:p("preset"),xyzBox:p("xyz-box"),xyz:p("xyz"),loadXyz:p("load-xyz"),charge:p("charge"),multiplicity:p("multiplicity"),method:p("method"),basis:p("basis"),relax:p("relax"),resetGeometry:p("reset-geometry"),showCloud:p("show-cloud"),cloudOpacity:p("cloud-opacity"),cloudMin:p("cloud-min"),showOrbital:p("show-orbital"),iso:p("iso"),status:p("status"),readouts:p("readouts"),orbitalName:p("orbital-name"),diagram:document.getElementById("diagram"),saveImage:p("save-image"),playScf:p("play-scf"),scfFrame:p("scf-frame"),scfLabel:p("scf-label"),help:p("help"),intro:p("intro"),legendCloud:p("legend-cloud"),legendOrbital:p("legend-orbital"),legendOrbitalName:p("legend-orbital-name"),legendAtoms:p("legend-atoms"),legendHelp:p("legend-help"),legendTitle:p("legend-title"),presetBlurb:p("preset-blurb"),settingsSummary:p("settings-summary"),energy:p("energy"),labels:p("labels"),tip:p("tip"),tipText:p("tip-text"),tipClose:p("tip-close"),restartTips:p("restart-tips"),introMore:p("intro-more"),settings:p("settings")},a={atoms:[],initialAtoms:[],charge:0,multiplicity:1,basis:"6-31g",method:"hf",result:null,orbital:null,frame:null,relaxId:null,hover:null},S=new nt,R=new Ke(new Worker(new URL(""+new URL("engine.worker-BsSZgK8L.js",import.meta.url).href,import.meta.url),{type:"module"}));let z=null;const ve=new URLSearchParams(location.search).has("headless"),Se=Math.max(16,Math.min(160,Number(new URLSearchParams(location.search).get("res"))||112));let A=!0,$=null;function ze(){return{showCloud:s.showCloud.checked,showOrbital:s.showOrbital.checked&&a.orbital!==null,cloudOpacity:40*Number(s.cloudOpacity.value)**2,logMin:Number(s.cloudMin.value),logMax:1,iso:Number(s.iso.value),orbitalOpacity:.55}}function E(t,e="ok"){s.status.textContent=t,s.status.className=`status ${e==="ok"?"":e}`}function se(t={}){return{atoms:a.atoms.map(e=>({symbol:e.symbol,xyz:e.xyz})),charge:a.charge,multiplicity:a.multiplicity,basis:a.basis,method:a.method,...t}}function I(){if(!z)return;const{spheres:t,bonds:e}=vt(a.atoms,a.hover);z.setAtoms(t,e),A=!0}function ae(){const t=a.result;if(!z||!t)return;const e=a.frame===null?null:t.history[a.frame];try{z.updateVolume(ct(t,a.orbital,e??t,Se))}catch(n){E(n instanceof Error?n.message:String(n),"error")}A=!0}function Ce(){z&&A&&(A=!1,Bt(),z.render(S,ze()),Et()),requestAnimationFrame(Ce)}function Et(){const t=a.atoms;for(;s.labels.children.length<t.length;)s.labels.append(document.createElement("span"));for(;s.labels.children.length>t.length;)s.labels.lastElementChild.remove();const e=S.viewProjection,n=s.canvas.getBoundingClientRect();t.forEach((o,r)=>{const i=s.labels.children[r],[l,c,h]=o.xyz,u=e[3]*l+e[7]*c+e[11]*h+e[15];if(i.textContent=o.symbol,i.hidden=u<=0,u<=0)return;const d=(e[0]*l+e[4]*c+e[8]*h+e[12])/u,m=(e[1]*l+e[5]*c+e[9]*h+e[13])/u;i.style.transform=`translate(${n.left+(d+1)/2*n.width}px, ${n.top+(1-m)/2*n.height}px) translate(-50%, -50%)`})}function Bt(){const t=Math.min(window.devicePixelRatio||1,2),e=Math.max(1,Math.floor(s.canvas.clientWidth*t)),n=Math.max(1,Math.floor(s.canvas.clientHeight*t));(s.canvas.width!==e||s.canvas.height!==n)&&(s.canvas.width=e,s.canvas.height=n)}function Ae(t){return{spin:"alpha",index:t.nalpha-1}}function le(t,e=!1){const n=a.result&&a.result.nmo===t.nmo&&a.result.unrestricted===t.unrestricted;a.result=t,(!n||!a.orbital||a.orbital.index>=t.nmo)&&(a.orbital=Ae(t)),e||H(),a.frame=null,s.scfFrame.max=String(Math.max(0,t.history.length-1)),s.scfFrame.value=s.scfFrame.max,s.scfFrame.disabled=t.history.length===0,s.scfLabel.textContent=t.history.length?`step ${t.history.length} of ${t.history.length}`:"calculation steps",s.energy.textContent=`${t.energy.toFixed(5)} Ha`,s.readouts.replaceChildren(...ut(t).flatMap(o=>{const r=document.createElement("dt");r.textContent=o.label;const i=document.createElement("dd");return i.textContent=o.value,o.title&&(i.title=o.title),[r,i]})),ce(),ae()}function ce(){const t=a.result;if(!t)return;const e=s.showOrbital.checked&&a.orbital!==null;if(bt(s.diagram,t,e?a.orbital:null,n=>{a.orbital=n,s.showOrbital.checked=!0,ce(),ae(),V("orbital")}),s.orbitalName.textContent="",e&&a.orbital){const n=Ee(t,a.orbital.spin)[a.orbital.index],o=`${n.label}${t.unrestricted?a.orbital.spin==="alpha"?" α":" β":""}`;s.orbitalName.textContent=`· ${o}`,s.legendOrbitalName.textContent=`${o} orbital`,s.legendOrbital.title=`${o}: ${pt(n.label)}. Blue where its wave is positive, orange where it is negative; where they meet, these electrons are never found.`}ue()}function ue(){s.legendTitle.textContent=s.preset.value==="custom"?"Your molecule":re.find(t=>t.id===s.preset.value)?.name??"",s.legendCloud.hidden=!s.showCloud.checked,s.legendOrbital.hidden=!(s.showOrbital.checked&&a.orbital),s.legendAtoms.replaceChildren(...ft(a.atoms).map(t=>{const e=document.createElement("span");e.className="legend-item",e.title=t.name;const n=document.createElement("span");return n.className="swatch",n.style.background=t.color,e.append(n,t.symbol),e}))}function Le(t){return`Solved in ${Math.round(t.millis)} ms (${t.iterations} steps)`}async function F(t={}){if(a.atoms.length!==0){E("Solving…","busy");try{const e=await R.computeLatest(se({history:t.history??!0}));if(!e)return;le(e),E(e.converged?Le(e):`The calculation did not settle after ${e.iterations} steps`,e.converged?"ok":"error")}catch(e){E(e instanceof Error?e.message:String(e),"error")}}}function Ue(t,e,n){_(),a.atoms=t,a.initialAtoms=t.map(i=>({symbol:i.symbol,xyz:[...i.xyz]})),a.charge=e,a.multiplicity=n,a.result=null,a.orbital=null,s.charge.value=String(e),s.multiplicity.value=String(n),s.xyz.value=qe(t).split(`
`).slice(2).join(`
`),he(),R.reset();const{center:o,radius:r}=xt(t);S.fit(o,r),Z(),s.presetBlurb.textContent=Ye[s.preset.value]??"Your own molecule. Drag its atoms, relax it, and look at its orbitals like any other.",I(),ue(),F()}function Te(t){const e=We(t);for(const n of Array.from(s.basis.options))n.disabled=!e.bases.includes(n.value);e.bases.includes(a.basis)||(a.basis=e.bases[e.bases.length-1],s.basis.value=a.basis),Ue(Me(e.xyz),e.charge,e.multiplicity)}function he(){const t=a.multiplicity!==1;for(const e of Array.from(s.method.options))e.disabled=t&&e.value!=="hf";t&&a.method!=="hf"&&(a.method="hf",s.method.value="hf")}function _(){a.relaxId!==null&&R.cancelRelax(a.relaxId),a.relaxId=null,s.relax.textContent="Relax geometry"}function Pt(){if(a.relaxId!==null)return _();s.relax.textContent="Stop",E("Relaxing…","busy"),a.relaxId=R.relax(se({history:!1}),{onStep:(t,e)=>{a.atoms=a.atoms.map((n,o)=>({symbol:n.symbol,xyz:t.coordinates.slice(3*o,3*o+3)})),I(),le(t),E(`Relaxing: step ${e+1}, E = ${t.energy.toFixed(6)} Ha, max force ${t.maxForce?.toExponential(1)}`,"busy")},onDone:(t,e)=>{a.relaxId=null,V("relax"),s.relax.textContent="Relax geometry",E(t?`Relaxed in ${e} steps · ${a.result.energy.toFixed(6)} Ha`:`Stopped after ${e} steps without converging`,t?"ok":"error")},onError:t=>{a.relaxId=null,s.relax.textContent="Relax geometry",E(t,"error")}})}function H(){$!==null&&clearInterval($),$=null,s.playScf.textContent="▶ Watch it solve"}function ne(t){const e=a.result;if(!e||e.history.length===0)return;a.frame=t>=e.history.length-1?null:t,s.scfFrame.value=String(t);const n=e.history[t];s.scfLabel.textContent=`step ${t+1} of ${e.history.length} · E = ${n.energy.toFixed(6)} Ha · error ${n.error.toExponential(1)}`,ae()}async function Ot(){if($!==null)return H();_(),E("Recording the calculation…","busy");const t=await R.compute(se({history:!0,warmStart:!1})).catch(n=>(E(n.message,"error"),null));if(!t)return;le(t,!0),E(Le(t));let e=0;s.playScf.textContent="■ Stop",ne(0),$=window.setInterval(()=>{if(e++,e>=t.history.length)return H();ne(e)},450)}function Re(t){const e=s.canvas.getBoundingClientRect();return[(t.clientX-e.left)/e.width*2-1,1-(t.clientY-e.top)/e.height*2]}function Fe(t){const{origin:e,dir:n}=S.ray(...Re(t));let o=null,r=1/0;return a.atoms.forEach((i,l)=>{const c=tt(e,n,i.xyz,Oe(i.symbol)*1.1);c!==null&&c<r&&(r=c,o=l)}),o}let w=null;s.canvas.addEventListener("pointerdown",t=>{s.canvas.setPointerCapture(t.pointerId);const e=Fe(t);e!==null?(_(),H(),w={kind:"atom",index:e,planePoint:[...a.atoms[e].xyz],normal:S.basis().forward}):w={kind:"orbit",x:t.clientX,y:t.clientY},s.canvas.classList.add("dragging")});s.canvas.addEventListener("pointermove",t=>{if(!w){const c=Fe(t);s.canvas.classList.toggle("over-atom",c!==null),c!==a.hover&&(a.hover=c,I());return}if(w.kind==="orbit"){S.orbit(t.clientX-w.x,t.clientY-w.y),w.x=t.clientX,w.y=t.clientY,A=!0;return}const{origin:e,dir:n}=S.ray(...Re(t)),o=T(n,w.normal);if(Math.abs(o)<1e-6)return;const r=T(W(w.planePoint,e),w.normal)/o,i=Be(e,oe(n,r)),l=w.index;a.atoms=a.atoms.map((c,h)=>h===l?{symbol:c.symbol,xyz:i}:c),I(),F({history:!1})});const Ne=()=>{w?.kind==="atom"&&V("drag"),w=null,s.canvas.classList.remove("dragging")};s.canvas.addEventListener("pointerup",Ne);s.canvas.addEventListener("pointercancel",Ne);s.canvas.addEventListener("wheel",t=>{t.preventDefault(),S.zoom(Math.exp(t.deltaY*.001)),A=!0},{passive:!1});window.addEventListener("resize",()=>A=!0);for(const t of re)s.preset.add(new Option(t.name,t.id));s.preset.add(new Option("Custom (pasted XYZ)","custom"));s.preset.addEventListener("change",()=>{s.preset.value==="custom"?(s.settings.open=!0,s.xyzBox.open=!0):(Te(s.preset.value),V("preset"))});s.loadXyz.addEventListener("click",()=>{try{const t=Me(s.xyz.value);if(t.length===0)throw new Error('No atoms found. Use lines like "O 0 0 0" in Angstrom.');s.preset.value="custom";for(const n of Array.from(s.basis.options))n.disabled=!1;const e=we({atoms:t,charge:a.charge,multiplicity:1});Ue(t,a.charge,e%2===0?1:2)}catch(t){E(t instanceof Error?t.message:String(t),"error")}});s.charge.addEventListener("change",()=>{a.charge=Math.round(Number(s.charge.value)),(we({atoms:a.atoms,charge:a.charge})+a.multiplicity-1)%2!==0&&(a.multiplicity=a.multiplicity===1?2:a.multiplicity-1),s.multiplicity.value=String(a.multiplicity),he(),F()});s.multiplicity.addEventListener("change",()=>{a.multiplicity=Math.max(1,Math.round(Number(s.multiplicity.value))),he(),F()});s.method.addEventListener("change",()=>{a.method=s.method.value,Z(),F()});s.basis.addEventListener("change",()=>{a.basis=s.basis.value,Z(),F()});s.relax.addEventListener("click",Pt);s.resetGeometry.addEventListener("click",()=>{_(),a.atoms=a.initialAtoms.map(t=>({symbol:t.symbol,xyz:[...t.xyz]})),I(),F()});for(const t of[s.cloudOpacity,s.cloudMin,s.iso])t.addEventListener("input",()=>A=!0);s.showCloud.addEventListener("input",()=>{ue(),A=!0});s.showOrbital.addEventListener("input",()=>{s.showOrbital.checked&&!a.orbital&&a.result&&(a.orbital=Ae(a.result)),s.showOrbital.checked&&V("orbital"),ce(),A=!0});s.playScf.addEventListener("click",Ot);s.saveImage.addEventListener("click",async()=>{const t=Math.min(window.devicePixelRatio||1,2),e=document.createElement("a");e.href=await Ge(Math.round(s.canvas.clientWidth*t),Math.round(s.canvas.clientHeight*t)),e.download=`${s.preset.value}-${a.method}-${a.basis}.png`,e.click()});s.scfFrame.addEventListener("input",()=>{H(),ne(Number(s.scfFrame.value))});const ke="quantum-molecules:intro-seen";function de(t){s.intro.open||s.intro.showModal(),s.tip.hidden=!0,t?(s.introMore.open=!0,document.getElementById(t)?.scrollIntoView({block:"start"})):s.intro.scrollTop=0}s.help.addEventListener("click",()=>de());s.legendHelp.addEventListener("click",()=>de("intro-colors"));s.intro.addEventListener("close",()=>{Y(ke),K()});s.intro.addEventListener("click",t=>{t.target===s.intro&&s.intro.close()});function Y(t,e=!0){try{e?localStorage.setItem(t,"1"):localStorage.removeItem(t)}catch{}}function xe(t){try{return localStorage.getItem(t)!==null}catch{return!1}}const X="quantum-molecules:tips-done",fe=new Set;let k=!1;function K(){const t=k?Mt(fe):null;t===null&&k&&(k=!1,Y(X)),s.tip.hidden=t===null||s.intro.open,t&&(s.tipText.innerHTML=t.html)}function V(t){fe.add(t),K()}s.tipClose.addEventListener("click",()=>{k=!1,Y(X),K()});s.restartTips.addEventListener("click",()=>{fe.clear(),k=!0,Y(X,!1),s.intro.close()});function Z(){s.settingsSummary.textContent=`(${s.method.selectedOptions[0]?.text}, ${s.basis.selectedOptions[0]?.text})`}async function St(){const t=!ve||new URLSearchParams(location.search).has("intro");k=t&&!xe(X),t&&!xe(ke)?de():K(),s.method.value=a.method,s.basis.value=a.basis,Z();try{z=await ie.create(ve?null:s.canvas),z.device.lost.then(e=>E(`The GPU device was lost: ${e.message}`,"error"))}catch(e){s.noGpu.hidden=!1,s.noGpu.textContent=`${e instanceof Error?e.message:String(e)} The calculations still run and the numbers on the right are real, but the 3D view needs a browser with WebGPU (recent Chrome, Edge, or Safari).`}s.preset.value="h2o",await R.ready,Te("h2o"),requestAnimationFrame(Ce)}async function Ge(t=s.canvas.clientWidth||800,e=s.canvas.clientHeight||600){if(!z)throw new Error("WebGPU is not available.");S.aspect=t/e;const n=await z.snapshot(S,ze(),t,e),o=document.createElement("canvas");return o.width=t,o.height=e,o.getContext("2d").putImageData(new ImageData(n,t,e),0,0),o.toDataURL("image/png")}Object.assign(window,{__qm:{state:a,engine:R,camera:S,snapshot:Ge,volumeResolution:Se,get renderer(){return z}}});St();
