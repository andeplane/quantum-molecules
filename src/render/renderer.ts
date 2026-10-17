import atomsWgsl from './atoms.wgsl?raw';
import raymarchWgsl from './raymarch.wgsl?raw';
import volumeWgsl from './volume.wgsl?raw';
import type { OrbitCamera } from './camera';
import { invert, type Vec3 } from './math';
import type { VolumeBox, VolumeInputs } from './volume-data';

export interface SphereInstance {
  center: Vec3;
  radius: number;
  color: Vec3;
  highlight: boolean;
}

export interface BondInstance {
  a: Vec3;
  b: Vec3;
  radius: number;
  colorA: Vec3;
  colorB: Vec3;
}

export interface ViewSettings {
  showCloud: boolean;
  showOrbital: boolean;
  /** Opacity per bohr of the densest part of the cloud. */
  cloudOpacity: number;
  /** log10 of the density shown as fully faint / fully bright. */
  logMin: number;
  logMax: number;
  iso: number;
  orbitalOpacity: number;
}

const DEPTH_FORMAT: GPUTextureFormat = 'depth24plus';
const RING_SEGMENTS = 18;

/** A unit cylinder as a triangle list of (cos, sin, along). */
function cylinderMesh(): Float32Array<ArrayBuffer> {
  const v: number[] = [];
  for (let i = 0; i < RING_SEGMENTS; i++) {
    const a0 = (i / RING_SEGMENTS) * 2 * Math.PI;
    const a1 = ((i + 1) / RING_SEGMENTS) * 2 * Math.PI;
    const p = (a: number, z: number) => v.push(Math.cos(a), Math.sin(a), z);
    p(a0, 0), p(a1, 0), p(a1, 1);
    p(a0, 0), p(a1, 1), p(a0, 1);
  }
  return new Float32Array(v);
}

export class Renderer {
  readonly device: GPUDevice;
  private context: GPUCanvasContext | null;
  readonly format: GPUTextureFormat;
  private depth: GPUTexture | null = null;
  private cameraBuffer: GPUBuffer;
  private rayBuffer: GPUBuffer;
  private paramsBuffer: GPUBuffer;
  private normsBuffer: GPUBuffer;
  private ringBuffer: GPUBuffer;
  private ringCount: number;
  private sphereBuffer: GPUBuffer | null = null;
  private sphereCount = 0;
  private bondBuffer: GPUBuffer | null = null;
  private bondCount = 0;
  private spherePipeline: GPURenderPipeline;
  private bondPipeline: GPURenderPipeline;
  private rayPipeline: GPURenderPipeline;
  private volumePipeline: GPUComputePipeline;
  private atomsBindGroup: GPUBindGroup;
  private sampler: GPUSampler;
  private volume: GPUTexture | null = null;
  private box: VolumeBox | null = null;
  private frame = 0;
  adapter: GPUAdapter | null = null;
  background: GPUColor = { r: 0.035, g: 0.045, b: 0.07, a: 1 };

  /** With `canvas` null the renderer only draws offscreen (see `snapshot`), which headless browsers need. */
  static async create(canvas: HTMLCanvasElement | null): Promise<Renderer> {
    if (!navigator.gpu) throw new Error('WebGPU is not available in this browser.');
    const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
    if (!adapter) throw new Error('No WebGPU adapter was found.');
    const device = await adapter.requestDevice();
    const context = canvas ? canvas.getContext('webgpu') : null;
    if (canvas && !context) throw new Error('Could not create a WebGPU canvas context.');
    const r = new Renderer(device, context, navigator.gpu.getPreferredCanvasFormat());
    // Chrome can lose the device if the adapter is garbage collected.
    r.adapter = adapter;
    return r;
  }

  constructor(device: GPUDevice, context: GPUCanvasContext | null, format: GPUTextureFormat) {
    this.device = device;
    this.context = context;
    this.format = format;
    context?.configure({ device, format, alphaMode: 'opaque' });

    const uniform = GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST;
    this.cameraBuffer = device.createBuffer({ size: 128, usage: uniform });
    this.rayBuffer = device.createBuffer({ size: 144, usage: uniform });
    this.paramsBuffer = device.createBuffer({ size: 48, usage: uniform });
    this.normsBuffer = device.createBuffer({ size: 40, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    const ring = cylinderMesh();
    this.ringCount = ring.length / 3;
    this.ringBuffer = device.createBuffer({ size: ring.byteLength, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
    device.queue.writeBuffer(this.ringBuffer, 0, ring);
    this.sampler = device.createSampler({ magFilter: 'linear', minFilter: 'linear' });

    const atomsModule = device.createShaderModule({ code: atomsWgsl, label: 'atoms' });
    const atomsLayout = device.createPipelineLayout({ bindGroupLayouts: [this.cameraLayout()] });
    const depthStencil: GPUDepthStencilState = { format: DEPTH_FORMAT, depthWriteEnabled: true, depthCompare: 'less' };
    const f32x3 = (shaderLocation: number, offset: number): GPUVertexAttribute => ({ shaderLocation, offset, format: 'float32x3' });
    const f32x1 = (shaderLocation: number, offset: number): GPUVertexAttribute => ({ shaderLocation, offset, format: 'float32' });
    this.spherePipeline = device.createRenderPipeline({
      label: 'spheres',
      layout: atomsLayout,
      vertex: {
        module: atomsModule,
        entryPoint: 'sphereVs',
        buffers: [{ arrayStride: 32, stepMode: 'instance', attributes: [f32x3(0, 0), f32x1(1, 12), f32x3(2, 16), f32x1(3, 28)] }],
      },
      fragment: { module: atomsModule, entryPoint: 'sphereFs', targets: [{ format }] },
      primitive: { topology: 'triangle-strip' },
      depthStencil,
    });
    this.bondPipeline = device.createRenderPipeline({
      label: 'bonds',
      layout: atomsLayout,
      vertex: {
        module: atomsModule,
        entryPoint: 'bondVs',
        buffers: [
          { arrayStride: 52, stepMode: 'instance', attributes: [f32x3(0, 0), f32x1(1, 12), f32x3(2, 16), f32x3(3, 28), f32x3(4, 40)] },
          { arrayStride: 12, stepMode: 'vertex', attributes: [f32x3(5, 0)] },
        ],
      },
      fragment: { module: atomsModule, entryPoint: 'bondFs', targets: [{ format }] },
      primitive: { topology: 'triangle-list', cullMode: 'back' },
      depthStencil,
    });
    this.atomsBindGroup = device.createBindGroup({
      layout: this.cameraLayout(),
      entries: [{ binding: 0, resource: { buffer: this.cameraBuffer } }],
    });

    const rayModule = device.createShaderModule({ code: raymarchWgsl, label: 'raymarch' });
    this.rayPipeline = device.createRenderPipeline({
      label: 'raymarch',
      layout: 'auto',
      vertex: { module: rayModule, entryPoint: 'vs' },
      fragment: {
        module: rayModule,
        entryPoint: 'fs',
        targets: [
          {
            format,
            blend: {
              color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' },
              alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' },
            },
          },
        ],
      },
    });

    this.volumePipeline = device.createComputePipeline({
      label: 'volume',
      layout: 'auto',
      compute: { module: device.createShaderModule({ code: volumeWgsl, label: 'volume' }), entryPoint: 'main' },
    });
  }

  private cameraLayoutCache: GPUBindGroupLayout | null = null;
  private cameraLayout(): GPUBindGroupLayout {
    this.cameraLayoutCache ??= this.device.createBindGroupLayout({
      entries: [{ binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } }],
    });
    return this.cameraLayoutCache;
  }

  resize(width: number, height: number) {
    if (this.depth && this.depth.width === width && this.depth.height === height) return;
    this.depth?.destroy();
    this.depth = this.device.createTexture({
      size: [width, height],
      format: DEPTH_FORMAT,
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING,
    });
  }

  setAtoms(spheres: SphereInstance[], bonds: BondInstance[]) {
    const s = new Float32Array(Math.max(1, spheres.length) * 8);
    spheres.forEach((a, i) => s.set([...a.center, a.radius, ...a.color, a.highlight ? 1 : 0], 8 * i));
    if (!this.sphereBuffer || this.sphereBuffer.size < s.byteLength) {
      this.sphereBuffer?.destroy();
      this.sphereBuffer = this.device.createBuffer({ size: s.byteLength, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
    }
    this.device.queue.writeBuffer(this.sphereBuffer, 0, s);
    this.sphereCount = spheres.length;

    const b = new Float32Array(Math.max(1, bonds.length) * 13);
    bonds.forEach((d, i) => b.set([...d.a, d.radius, ...d.b, ...d.colorA, ...d.colorB], 13 * i));
    if (!this.bondBuffer || this.bondBuffer.size < b.byteLength) {
      this.bondBuffer?.destroy();
      this.bondBuffer = this.device.createBuffer({ size: b.byteLength, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
    }
    this.device.queue.writeBuffer(this.bondBuffer, 0, b);
    this.bondCount = bonds.length;
  }

  /** Evaluate the density and the selected orbital on the GPU grid. */
  updateVolume(v: VolumeInputs) {
    const { device } = this;
    const [nx, ny, nz] = v.box.dims;
    if (!this.volume || this.volume.width !== nx || this.volume.height !== ny || this.volume.depthOrArrayLayers !== nz) {
      this.volume?.destroy();
      this.volume = device.createTexture({
        size: [nx, ny, nz],
        dimension: '3d',
        format: 'rgba16float',
        usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC,
      });
    }
    this.box = v.box;
    const params = new ArrayBuffer(48);
    const f = new Float32Array(params);
    const u = new Uint32Array(params);
    f.set(v.box.min, 0);
    u[3] = v.nbf;
    f.set(v.box.max, 4);
    u[7] = v.nshell;
    u.set(v.box.dims, 8);
    u[11] = v.nocc;
    device.queue.writeBuffer(this.paramsBuffer, 0, params);
    device.queue.writeBuffer(this.normsBuffer, 0, v.norms);
    const storage = (data: ArrayBuffer | Float32Array<ArrayBuffer>) => {
      const buf = device.createBuffer({ size: Math.max(16, data.byteLength), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
      device.queue.writeBuffer(buf, 0, data);
      return buf;
    };
    const buffers = [storage(v.shells), storage(v.prims), storage(v.coefs), storage(v.weights.length ? v.weights : new Float32Array(4))];
    const bindGroup = device.createBindGroup({
      layout: this.volumePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.paramsBuffer } },
        { binding: 1, resource: { buffer: buffers[0] } },
        { binding: 2, resource: { buffer: buffers[1] } },
        { binding: 3, resource: { buffer: this.normsBuffer } },
        { binding: 4, resource: { buffer: buffers[2] } },
        { binding: 5, resource: { buffer: buffers[3] } },
        { binding: 6, resource: this.volume.createView() },
      ],
    });
    const enc = device.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(this.volumePipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(nx / 4), Math.ceil(ny / 4), Math.ceil(nz / 4));
    pass.end();
    device.queue.submit([enc.finish()]);
    for (const b of buffers) b.destroy();
  }

  /** Draw to the canvas, or to `target` (a texture in `this.format`). */
  render(camera: OrbitCamera, view: ViewSettings, target?: GPUTexture) {
    const { device } = this;
    const tex = target ?? this.context?.getCurrentTexture();
    if (!tex) return;
    this.resize(tex.width, tex.height);
    camera.aspect = tex.width / tex.height;
    const cam = new Float32Array(32);
    cam.set(camera.view, 0);
    cam.set(camera.projection, 16);
    device.queue.writeBuffer(this.cameraBuffer, 0, cam);

    const enc = device.createCommandEncoder();
    const atomsPass = enc.beginRenderPass({
      colorAttachments: [{ view: tex.createView(), clearValue: this.background, loadOp: 'clear', storeOp: 'store' }],
      depthStencilAttachment: { view: this.depth!.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
    });
    if (this.bondCount > 0 && this.bondBuffer) {
      atomsPass.setPipeline(this.bondPipeline);
      atomsPass.setBindGroup(0, this.atomsBindGroup);
      atomsPass.setVertexBuffer(0, this.bondBuffer);
      atomsPass.setVertexBuffer(1, this.ringBuffer);
      atomsPass.draw(this.ringCount, this.bondCount);
    }
    if (this.sphereCount > 0 && this.sphereBuffer) {
      atomsPass.setPipeline(this.spherePipeline);
      atomsPass.setBindGroup(0, this.atomsBindGroup);
      atomsPass.setVertexBuffer(0, this.sphereBuffer);
      atomsPass.draw(4, this.sphereCount);
    }
    atomsPass.end();

    if (this.volume && this.box && (view.showCloud || view.showOrbital)) {
      const box = this.box;
      const u = new Float32Array(36);
      u.set(invert(camera.viewProjection), 0);
      u.set(camera.eye, 16);
      u[19] = view.cloudOpacity;
      u.set(box.min, 20);
      u[23] = view.logMin;
      u.set(box.max, 24);
      u[27] = view.logMax;
      u.set(
        box.max.map((m, k) => (m - box.min[k]) / Math.max(1, box.dims[k] - 1)),
        28,
      );
      u[31] = view.iso;
      u[32] = view.showCloud ? 1 : 0;
      u[33] = view.showOrbital ? 1 : 0;
      u[34] = view.orbitalOpacity;
      u[35] = this.frame++ % 64;
      device.queue.writeBuffer(this.rayBuffer, 0, u);
      const bindGroup = device.createBindGroup({
        layout: this.rayPipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: this.rayBuffer } },
          { binding: 1, resource: this.volume.createView() },
          { binding: 2, resource: this.sampler },
          { binding: 3, resource: this.depth!.createView() },
        ],
      });
      const rayPass = enc.beginRenderPass({
        colorAttachments: [{ view: tex.createView(), loadOp: 'load', storeOp: 'store' }],
      });
      rayPass.setPipeline(this.rayPipeline);
      rayPass.setBindGroup(0, bindGroup);
      rayPass.draw(3);
      rayPass.end();
    }
    device.queue.submit([enc.finish()]);
  }

  /** Render offscreen and return RGBA pixels, for saving images and for headless tests. */
  async snapshot(camera: OrbitCamera, view: ViewSettings, width: number, height: number): Promise<Uint8ClampedArray<ArrayBuffer>> {
    const texture = this.device.createTexture({
      size: [width, height],
      format: this.format,
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
    });
    this.render(camera, view, texture);
    const bytesPerRow = Math.ceil((width * 4) / 256) * 256;
    const buf = this.device.createBuffer({ size: bytesPerRow * height, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
    const enc = this.device.createCommandEncoder();
    enc.copyTextureToBuffer({ texture }, { buffer: buf, bytesPerRow }, [width, height]);
    this.device.queue.submit([enc.finish()]);
    await buf.mapAsync(GPUMapMode.READ);
    const src = new Uint8Array(buf.getMappedRange());
    const out = new Uint8ClampedArray(width * height * 4);
    const bgra = this.format.startsWith('bgra');
    for (let y = 0; y < height; y++)
      for (let x = 0; x < width; x++) {
        const i = y * bytesPerRow + 4 * x;
        const o = 4 * (y * width + x);
        out[o] = src[i + (bgra ? 2 : 0)];
        out[o + 1] = src[i + 1];
        out[o + 2] = src[i + (bgra ? 0 : 2)];
        out[o + 3] = 255;
      }
    buf.unmap();
    buf.destroy();
    texture.destroy();
    return out;
  }

  /** Read the volume back as [rho, psi] pairs in x-fastest order (for tests and debugging). */
  async readVolume(): Promise<{ dims: [number, number, number]; data: Float32Array } | null> {
    if (!this.volume || !this.box) return null;
    const [nx, ny, nz] = this.box.dims;
    const bytesPerRow = Math.ceil((nx * 8) / 256) * 256;
    const buf = this.device.createBuffer({ size: bytesPerRow * ny * nz, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
    const enc = this.device.createCommandEncoder();
    enc.copyTextureToBuffer({ texture: this.volume }, { buffer: buf, bytesPerRow, rowsPerImage: ny }, [nx, ny, nz]);
    this.device.queue.submit([enc.finish()]);
    await buf.mapAsync(GPUMapMode.READ);
    const half = new Uint16Array(buf.getMappedRange());
    const data = new Float32Array(nx * ny * nz * 2);
    for (let z = 0; z < nz; z++)
      for (let y = 0; y < ny; y++)
        for (let x = 0; x < nx; x++) {
          const src = ((z * ny + y) * bytesPerRow) / 2 + x * 4;
          const dst = 2 * ((z * ny + y) * nx + x);
          data[dst] = halfToFloat(half[src]);
          data[dst + 1] = halfToFloat(half[src + 1]);
        }
    buf.unmap();
    buf.destroy();
    return { dims: [nx, ny, nz], data };
  }
}

export function halfToFloat(h: number): number {
  const s = h & 0x8000 ? -1 : 1;
  const e = (h >> 10) & 0x1f;
  const m = h & 0x3ff;
  if (e === 0) return s * m * 2 ** -24;
  if (e === 31) return m ? NaN : s * Infinity;
  return s * (1 + m / 1024) * 2 ** (e - 15);
}
