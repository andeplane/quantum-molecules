import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { parseXyz } from '../src/chem/molecule';
import { BLURBS, PRESETS, preset } from '../src/chem/presets';
import { TIPS, currentTip } from '../src/ui/tips';
import { CLOUD_RAMP, LOBE_NEGATIVE, LOBE_POSITIVE, css, describeOrbital, legendElements } from '../src/ui/legend';

const shader = readFileSync(new URL('../src/render/raymarch.wgsl', import.meta.url), 'utf8');
const styles = readFileSync(new URL('../src/style.css', import.meta.url), 'utf8');
const vec3 = (src: string) => src.split(',').map(Number);

describe('legend', () => {
  it('lists each element once, heaviest first', () => {
    const water = legendElements(parseXyz(preset('h2o').xyz));
    expect(water).toEqual([
      { symbol: 'O', name: 'oxygen', color: 'rgb(255, 33, 33)' },
      { symbol: 'H', name: 'hydrogen', color: 'rgb(242, 242, 242)' },
    ]);
    expect(legendElements(parseXyz(preset('ch2o').xyz)).map((e) => e.symbol)).toEqual(['O', 'C', 'H']);
    expect(legendElements(parseXyz('h 0 0 0\nH 0 0 0.74'))).toHaveLength(1);
    expect(legendElements([])).toEqual([]);
  });

  it('describes orbitals in plain words', () => {
    expect(describeOrbital('HOMO')).toContain('highest-energy orbital that holds electrons');
    expect(describeOrbital('LUMO')).toContain('lowest-energy empty orbital');
    expect(describeOrbital('HOMO−2')).toContain('occupied');
    expect(describeOrbital('LUMO+1')).toContain('empty');
  });

  it('matches the colours the shader draws', () => {
    const lobes = /select\(vec3f\(([^)]*)\), vec3f\(([^)]*)\), level > 0\.0\)/.exec(shader);
    expect(lobes).not.toBeNull();
    expect(vec3(lobes![1])).toEqual(LOBE_NEGATIVE);
    expect(vec3(lobes![2])).toEqual(LOBE_POSITIVE);
    const ramp = [...shader.matchAll(/let c[012] = vec3f\(([^)]*)\);/g)].map((m) => vec3(m[1]));
    expect(ramp).toEqual(CLOUD_RAMP);
  });

  it('uses the same colours for the swatches', () => {
    expect(styles).toContain(`.lobe-pos {\n  background: ${css(LOBE_POSITIVE)};`);
    expect(styles).toContain(`.lobe-neg {\n  background: ${css(LOBE_NEGATIVE)};`);
    expect(styles).toContain(`linear-gradient(90deg, ${CLOUD_RAMP.map(css).join(', ')})`);
  });

  it('formats CSS colours', () => {
    expect(css([1, 0.5, 0])).toBe('rgb(255, 128, 0)');
  });
});

describe('newcomer text', () => {
  it('has a blurb for every preset', () => {
    for (const p of PRESETS) expect(BLURBS[p.id], p.id).toMatch(/\.$/);
  });

  it('walks through the tips in order, skipping what was already done', () => {
    expect(currentTip(new Set())?.until).toBe('drag');
    expect(currentTip(new Set(['drag']))?.until).toBe('relax');
    expect(currentTip(new Set(['drag', 'orbital']))?.until).toBe('relax');
    expect(currentTip(new Set(['drag', 'relax', 'orbital']))?.until).toBe('preset');
    expect(currentTip(new Set(TIPS.map((t) => t.until)))).toBeNull();
  });
});
