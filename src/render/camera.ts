import { add, invert, lookAt, multiply, normalize, perspective, scale, sub, transformPoint, type Mat4, type Vec3 } from './math';

/** Orbit camera around a target, in bohr. */
export class OrbitCamera {
  target: Vec3 = [0, 0, 0];
  distance = 15;
  yaw = 0.6;
  pitch = 0.35;
  fovy = (40 * Math.PI) / 180;
  aspect = 1;

  get eye(): Vec3 {
    const cp = Math.cos(this.pitch);
    return add(this.target, scale([cp * Math.sin(this.yaw), Math.sin(this.pitch), cp * Math.cos(this.yaw)], this.distance));
  }

  get view(): Mat4 {
    return lookAt(this.eye, this.target, [0, 1, 0]);
  }

  get projection(): Mat4 {
    return perspective(this.fovy, this.aspect, Math.max(0.05, this.distance * 0.02), this.distance * 10 + 50);
  }

  get viewProjection(): Mat4 {
    return multiply(this.projection, this.view);
  }

  orbit(dx: number, dy: number) {
    this.yaw -= dx * 0.008;
    this.pitch = Math.max(-1.5, Math.min(1.5, this.pitch + dy * 0.008));
  }

  zoom(factor: number) {
    this.distance = Math.max(2, Math.min(200, this.distance * factor));
  }

  /** Frame a bounding sphere. */
  fit(center: Vec3, radius: number) {
    this.target = center;
    this.distance = Math.max(6, (radius + 2.5) / Math.sin(this.fovy / 2));
  }

  /** World-space ray through a point in normalized device coordinates. */
  ray(ndcX: number, ndcY: number): { origin: Vec3; dir: Vec3 } {
    const inv = invert(this.viewProjection);
    const near = transformPoint(inv, [ndcX, ndcY, 0]);
    const far = transformPoint(inv, [ndcX, ndcY, 1]);
    return { origin: near, dir: normalize(sub(far, near)) };
  }

  /** Camera right and up vectors, for dragging in the view plane. */
  basis(): { right: Vec3; up: Vec3; forward: Vec3 } {
    const v = this.view;
    return { right: [v[0], v[4], v[8]], up: [v[1], v[5], v[9]], forward: [-v[2], -v[6], -v[10]] };
  }
}
