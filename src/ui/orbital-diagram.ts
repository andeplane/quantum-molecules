import { EV_PER_HARTREE } from '../chem/elements';
import { orbitalLevels } from '../chem/orbitals';
import type { ComputeOutput } from '../engine/types';
import type { OrbitalSelection } from '../render/volume-data';

export interface DiagramLevel {
  spin: 'alpha' | 'beta';
  index: number;
  energy: number;
  occupied: boolean;
  label: string;
  /** Vertical position in [0, 1], 0 at the top. */
  y: number;
}

/** Core levels this far (Hartree) below the HOMO are left out so the valence levels stay readable. */
export const CORE_WINDOW = 1.6;

/** Levels in a window around the frontier orbitals, with positions for drawing. */
export function diagramLevels(res: ComputeOutput, below = 8, above = 6): DiagramLevel[] {
  const spins: ('alpha' | 'beta')[] = res.unrestricted ? ['alpha', 'beta'] : ['alpha'];
  const levels = spins.flatMap((spin) => {
    const nocc = spin === 'alpha' ? res.nalpha : res.nbeta;
    const all = orbitalLevels(res, spin);
    const homo = all[nocc - 1]?.energy ?? all[0].energy;
    return all
      .filter((l) => l.index >= nocc - below && l.index <= nocc + above - 1 && l.energy > homo - CORE_WINDOW)
      .map((l) => ({ ...l, spin }));
  });
  const es = levels.map((l) => l.energy);
  const lo = Math.min(...es);
  const hi = Math.max(...es);
  const span = hi - lo || 1;
  return levels.map((l) => ({ ...l, y: 0.04 + 0.92 * ((hi - l.energy) / span) }));
}

const SVG_NS = 'http://www.w3.org/2000/svg';

export function renderDiagram(
  svg: SVGSVGElement,
  res: ComputeOutput,
  selected: OrbitalSelection | null,
  onSelect: (s: OrbitalSelection) => void,
) {
  svg.replaceChildren();
  const width = 240;
  const height = 300;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  const levels = diagramLevels(res);
  const columns = res.unrestricted ? 2 : 1;
  const colWidth = 70;
  const x0 = (s: 'alpha' | 'beta') => (columns === 1 ? 70 : s === 'alpha' ? 40 : 125);
  // Degenerate levels sit side by side.
  const groups = new Map<string, DiagramLevel[]>();
  for (const l of levels) {
    const key = `${l.spin}:${Math.round(l.energy * EV_PER_HARTREE * 50)}`;
    groups.set(key, [...(groups.get(key) ?? []), l]);
  }
  for (const group of groups.values()) {
    const w = (colWidth - 6 * (group.length - 1)) / group.length;
    group.forEach((l, k) => {
      const x = x0(l.spin) + k * (w + 6);
      const y = l.y * (height - 20) + 10;
      const isSel = selected?.spin === l.spin && selected.index === l.index;
      const g = document.createElementNS(SVG_NS, 'g');
      g.classList.add('level');
      if (l.occupied) g.classList.add('occupied');
      if (isSel) g.classList.add('selected');
      const hit = document.createElementNS(SVG_NS, 'rect');
      hit.setAttribute('x', String(x - 2));
      hit.setAttribute('y', String(y - 7));
      hit.setAttribute('width', String(w + 4));
      hit.setAttribute('height', '14');
      hit.setAttribute('class', 'hit');
      const line = document.createElementNS(SVG_NS, 'line');
      line.setAttribute('x1', String(x));
      line.setAttribute('x2', String(x + w));
      line.setAttribute('y1', String(y));
      line.setAttribute('y2', String(y));
      const title = document.createElementNS(SVG_NS, 'title');
      title.textContent = `${l.label}${columns === 2 ? ` (${l.spin === 'alpha' ? 'α' : 'β'})` : ''}: ${(l.energy * EV_PER_HARTREE).toFixed(2)} eV`;
      g.append(hit, line, title);
      if (l.occupied) {
        const n = res.unrestricted ? 1 : 2;
        for (let e = 0; e < n; e++) {
          const arrow = document.createElementNS(SVG_NS, 'text');
          arrow.setAttribute('x', String(x + w / 2 + (n === 2 ? (e === 0 ? -5 : 5) : 0)));
          arrow.setAttribute('y', String(y + 4));
          arrow.setAttribute('class', 'electron');
          arrow.textContent = (l.spin === 'beta') !== (e === 1) ? '↓' : '↑';
          g.append(arrow);
        }
      }
      if (k === group.length - 1 && (l.label === 'HOMO' || l.label === 'LUMO' || isSel)) {
        const text = document.createElementNS(SVG_NS, 'text');
        text.setAttribute('x', String(x + w + 6));
        text.setAttribute('y', String(y + 4));
        text.setAttribute('class', 'label');
        text.textContent = l.label;
        g.append(text);
      }
      g.addEventListener('click', () => onSelect({ spin: l.spin, index: l.index }));
      svg.append(g);
    });
  }
  if (columns === 2) {
    for (const [s, label] of [['alpha', 'α'], ['beta', 'β']] as const) {
      const t = document.createElementNS(SVG_NS, 'text');
      t.setAttribute('x', String(x0(s) + colWidth / 2));
      t.setAttribute('y', String(height - 2));
      t.setAttribute('class', 'spin-label');
      t.textContent = label;
      svg.append(t);
    }
  }
}
