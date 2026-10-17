import { BOHR_PER_ANGSTROM, element } from './elements';

export type Vec3 = [number, number, number];

export interface Atom {
  symbol: string;
  /** Position in bohr. */
  xyz: Vec3;
}

export interface Molecule {
  atoms: Atom[];
  charge: number;
  multiplicity: number;
}

export function nuclearCharge(atom: Atom): number {
  return element(atom.symbol).z;
}

export function electronCount(mol: Molecule): number {
  return mol.atoms.reduce((n, a) => n + nuclearCharge(a), 0) - mol.charge;
}

export function nuclearRepulsion(atoms: Atom[]): number {
  let e = 0;
  for (let i = 0; i < atoms.length; i++) {
    for (let j = 0; j < i; j++) {
      const [x, y, z] = atoms[i].xyz;
      const [u, v, w] = atoms[j].xyz;
      e += (nuclearCharge(atoms[i]) * nuclearCharge(atoms[j])) / Math.hypot(x - u, y - v, z - w);
    }
  }
  return e;
}

export function nuclearRepulsionGradient(atoms: Atom[]): Float64Array {
  const g = new Float64Array(atoms.length * 3);
  for (let i = 0; i < atoms.length; i++) {
    for (let j = 0; j < i; j++) {
      const d = [0, 1, 2].map((k) => atoms[i].xyz[k] - atoms[j].xyz[k]);
      const r = Math.hypot(d[0], d[1], d[2]);
      const f = (-nuclearCharge(atoms[i]) * nuclearCharge(atoms[j])) / (r * r * r);
      for (let k = 0; k < 3; k++) {
        g[3 * i + k] += f * d[k];
        g[3 * j + k] -= f * d[k];
      }
    }
  }
  return g;
}

/** Parse an XYZ block (Angstrom). The first two lines may be the usual count and comment lines. */
export function parseXyz(text: string): Atom[] {
  const atoms: Atom[] = [];
  for (const line of text.split(/\r?\n/)) {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) continue;
    const [sym, ...rest] = parts;
    const coords = rest.slice(0, 3).map(Number);
    if (!/^[A-Za-z]{1,2}$/.test(sym) || coords.some((c) => !Number.isFinite(c))) continue;
    const symbol = sym[0].toUpperCase() + sym.slice(1).toLowerCase();
    element(symbol);
    atoms.push({ symbol, xyz: coords.map((c) => c * BOHR_PER_ANGSTROM) as Vec3 });
  }
  return atoms;
}

export function toXyz(atoms: Atom[], comment = ''): string {
  const lines = atoms.map(
    (a) => `${a.symbol.padEnd(2)} ${a.xyz.map((c) => (c / BOHR_PER_ANGSTROM).toFixed(6).padStart(12)).join(' ')}`,
  );
  return `${atoms.length}\n${comment}\n${lines.join('\n')}\n`;
}
