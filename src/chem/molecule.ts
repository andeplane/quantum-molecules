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

/** Bond if closer than the scaled sum of covalent radii. Pairs of atom indices. */
export function bonds(atoms: Atom[], scale = 1.25): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < atoms.length; i++) {
    for (let j = 0; j < i; j++) {
      const [x, y, z] = atoms[i].xyz;
      const [u, v, w] = atoms[j].xyz;
      const r = Math.hypot(x - u, y - v, z - w) / BOHR_PER_ANGSTROM;
      if (r < scale * (element(atoms[i].symbol).covalentRadius + element(atoms[j].symbol).covalentRadius)) out.push([j, i]);
    }
  }
  return out;
}

/** Shift atoms so their geometric centre is at the origin. */
export function centered(atoms: Atom[]): Atom[] {
  const c = [0, 1, 2].map((k) => atoms.reduce((s, a) => s + a.xyz[k], 0) / Math.max(atoms.length, 1));
  return atoms.map((a) => ({ symbol: a.symbol, xyz: [a.xyz[0] - c[0], a.xyz[1] - c[1], a.xyz[2] - c[2]] }));
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
