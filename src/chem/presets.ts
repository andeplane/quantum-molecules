import type { BasisName } from '../engine/types';
import { BOHR_PER_ANGSTROM } from './elements';
import { CRYSTALS, buildCrystal } from './crystals';

export type PresetGroup = 'small' | 'large' | 'crystal';

export const GROUP_NAMES: Record<PresetGroup, string> = {
  small: 'Small molecules',
  large: 'Larger molecules',
  crystal: 'Pieces of crystals',
};

export interface Preset {
  id: string;
  name: string;
  group: PresetGroup;
  /** XYZ body in Angstrom (symbol x y z per line). */
  xyz: string;
  charge: number;
  multiplicity: number;
  /** Basis sets that stay interactive for this molecule. */
  bases: BasisName[];
  /** Key into `CRYSTALS` for a cluster cut from a crystal. */
  crystal?: string;
}

const ALL: BasisName[] = ['sto-3g', '6-31g', '6-31g*'];

function benzene(): string {
  const lines: string[] = [];
  for (let k = 0; k < 6; k++) {
    const t = (k * Math.PI) / 3;
    lines.push(`C ${(1.397 * Math.cos(t)).toFixed(4)} ${(1.397 * Math.sin(t)).toFixed(4)} 0`);
  }
  for (let k = 0; k < 6; k++) {
    const t = (k * Math.PI) / 3;
    lines.push(`H ${(2.481 * Math.cos(t)).toFixed(4)} ${(2.481 * Math.sin(t)).toFixed(4)} 0`);
  }
  return lines.join('\n');
}

/** The computed cluster of a crystal as an XYZ body (Angstrom). */
function crystalXyz(id: string): string {
  return buildCrystal(CRYSTALS[id])
    .cluster.map((a) => `${a.symbol} ${a.xyz.map((c) => (c / BOHR_PER_ANGSTROM).toFixed(4)).join(' ')}`)
    .join('\n');
}

const STO: BasisName[] = ['sto-3g'];

const GLYCINE =
  'N 1.2823 -0.5147 -0.5204\nC 0.5109 0.5487 0.1775\nC -0.9994 0.3600 -0.0969\nO -1.7660 1.1859 -0.5529\nO -1.4161 -0.9179 0.2628\nH 2.2759 -0.3804 -0.2705\nH 1.0177 -1.4145 -0.0875\nH 0.6645 0.5519 1.2627\nH 0.8096 1.5223 -0.2117\nH -2.3793 -0.9414 0.0368';
const PYRIDINE =
  'C 0.0000 -1.1889 -0.0000\nC 1.1902 -0.4804 0.0000\nC 1.1443 0.9069 -0.0000\nN -0.0000 1.6294 -0.0000\nC -1.1443 0.9069 0.0000\nC -1.1902 -0.4804 0.0000\nH 0.0000 -2.2723 -0.0000\nH 2.1430 -0.9925 0.0000\nH 2.0665 1.4821 0.0000\nH -2.0666 1.4820 0.0000\nH -2.1430 -0.9926 -0.0000';
const ASPIRIN =
  'C -3.4689 -0.5199 -0.0085\nC -2.0025 -0.4696 -0.4637\nO -1.5908 -0.3051 -1.5926\nO -1.1518 -0.6764 0.6362\nC 0.2190 -0.5834 0.3187\nC 0.9646 -1.7585 0.2673\nC 2.3291 -1.7078 0.0362\nC 2.9623 -0.4854 -0.1413\nC 2.2263 0.6857 -0.0824\nC 0.8535 0.6546 0.1471\nC 0.0563 1.9438 0.1941\nO -1.1555 2.0453 0.2931\nO 0.8687 3.0708 0.1013\nH -3.6626 -1.4583 0.5025\nH -4.1210 -0.4313 -0.8712\nH -3.6613 0.2977 0.6811\nH 0.4589 -2.7051 0.4067\nH 2.9011 -2.6268 -0.0081\nH 4.0277 -0.4466 -0.3281\nH 2.7107 1.6450 -0.2148\nH 0.2362 3.8312 0.1263';
const GLUCOSE =
  'O -3.3581 -0.4226 0.8698\nC -2.4381 0.1337 -0.0805\nC -1.0469 0.2517 0.6038\nO -0.5215 -1.0443 0.9703\nC 0.3704 -1.5787 -0.0096\nO -0.2886 -1.6403 -1.2813\nC 1.7078 -0.7394 -0.0736\nO 2.6939 -1.2430 0.8367\nC 1.4009 0.7770 0.1695\nO 2.2495 1.6487 -0.5896\nC -0.0595 1.0996 -0.2580\nO -0.3711 2.4826 -0.0486\nH -4.2037 -0.5061 0.3595\nH -2.3469 -0.4941 -0.9735\nH -2.7445 1.1386 -0.3973\nH -1.1910 0.7830 1.5517\nH 0.6257 -2.5902 0.3491\nH 0.2821 -2.2417 -1.8210\nH 2.1345 -0.8470 -1.0765\nH 2.2483 -1.2037 1.7211\nH 1.4974 0.9816 1.2444\nH 3.1598 1.4722 -0.2436\nH -0.1715 0.8305 -1.3158\nH 0.3712 2.9520 -0.5070';


/** Experimental equilibrium geometries (rounded). */
const SMALL: Omit<Preset, 'group'>[] = [
  { id: 'h2', name: 'Hydrogen H₂', xyz: 'H 0 0 -0.3705\nH 0 0 0.3705', charge: 0, multiplicity: 1, bases: ALL },
  {
    id: 'h2o',
    name: 'Water H₂O',
    xyz: 'O 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'nh3',
    name: 'Ammonia NH₃',
    xyz: 'N 0 0 0.1162\nH 0 0.9377 -0.2711\nH 0.8121 -0.4689 -0.2711\nH -0.8121 -0.4689 -0.2711',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'ch4',
    name: 'Methane CH₄',
    xyz: 'C 0 0 0\nH 0.6291 0.6291 0.6291\nH -0.6291 -0.6291 0.6291\nH -0.6291 0.6291 -0.6291\nH 0.6291 -0.6291 -0.6291',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  { id: 'hf', name: 'Hydrogen fluoride HF', xyz: 'F 0 0 0.0917\nH 0 0 -0.8253', charge: 0, multiplicity: 1, bases: ALL },
  { id: 'n2', name: 'Nitrogen N₂', xyz: 'N 0 0 -0.5488\nN 0 0 0.5488', charge: 0, multiplicity: 1, bases: ALL },
  { id: 'co2', name: 'Carbon dioxide CO₂', xyz: 'C 0 0 0\nO 0 0 1.1621\nO 0 0 -1.1621', charge: 0, multiplicity: 1, bases: ALL },
  {
    id: 'ch2o',
    name: 'Formaldehyde CH₂O',
    xyz: 'C 0 0 -0.5285\nO 0 0 0.6745\nH 0 0.9363 -1.1079\nH 0 -0.9363 -1.1079',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  {
    id: 'c2h4',
    name: 'Ethylene C₂H₄',
    xyz: 'C 0 0 0.6695\nC 0 0 -0.6695\nH 0 0.9288 1.2321\nH 0 -0.9288 1.2321\nH 0 0.9288 -1.2321\nH 0 -0.9288 -1.2321',
    charge: 0,
    multiplicity: 1,
    bases: ALL,
  },
  { id: 'c6h6', name: 'Benzene C₆H₆', xyz: benzene(), charge: 0, multiplicity: 1, bases: ['sto-3g'] },
  { id: 'oh', name: 'Hydroxyl radical OH', xyz: 'O 0 0 0\nH 0 0 0.9697', charge: 0, multiplicity: 2, bases: ALL },
  { id: 'o2', name: 'Oxygen O₂ (triplet)', xyz: 'O 0 0 -0.6038\nO 0 0 0.6038', charge: 0, multiplicity: 3, bases: ALL },
];

/** Geometries relaxed with this engine at RHF/STO-3G, starting from RDKit conformers. */
const LARGE: Omit<Preset, 'group'>[] = [
  { id: 'glycine', name: 'Glycine C₂H₅NO₂', xyz: GLYCINE, charge: 0, multiplicity: 1, bases: ALL },
  { id: 'pyridine', name: 'Pyridine C₅H₅N', xyz: PYRIDINE, charge: 0, multiplicity: 1, bases: ALL },
  { id: 'aspirin', name: 'Aspirin C₉H₈O₄', xyz: ASPIRIN, charge: 0, multiplicity: 1, bases: STO },
  { id: 'glucose', name: 'Glucose C₆H₁₂O₆', xyz: GLUCOSE, charge: 0, multiplicity: 1, bases: STO },
];

/** Clusters cut from the crystal lattice (see `crystals.ts`); atoms sit on their lattice sites. */
const CRYSTAL: Omit<Preset, 'group' | 'xyz'>[] = [
  { id: 'diamond', name: 'Diamond (C₁₀H₁₆ piece)', crystal: 'diamond', charge: 0, multiplicity: 1, bases: STO },
  { id: 'silicon', name: 'Silicon (Si₁₀H₁₆ piece)', crystal: 'silicon', charge: 0, multiplicity: 1, bases: STO },
  { id: 'graphite', name: 'Graphite (C₁₆H₁₀ flake)', crystal: 'graphite', charge: 0, multiplicity: 1, bases: STO },
  { id: 'nacl', name: 'Table salt NaCl (Na₄Cl₄ cube)', crystal: 'nacl', charge: 0, multiplicity: 1, bases: STO },
  { id: 'mgo', name: 'Magnesium oxide MgO (Mg₄O₄ cube)', crystal: 'mgo', charge: 0, multiplicity: 1, bases: STO },
];

export const PRESETS: Preset[] = [
  ...SMALL.map((p) => ({ ...p, group: 'small' as const })),
  ...LARGE.map((p) => ({ ...p, group: 'large' as const })),
  ...CRYSTAL.map((p) => ({ ...p, group: 'crystal' as const, xyz: crystalXyz(p.crystal!) })),
];

export function preset(id: string): Preset {
  const p = PRESETS.find((p) => p.id === id);
  if (!p) throw new Error(`Unknown preset ${id}`);
  return p;
}

/** One or two plain-language sentences about each preset, for people new to chemistry. */
export const BLURBS: Record<string, string> = {
  h2: 'The simplest molecule: two hydrogen atoms sharing their two electrons in one bond.',
  h2o: 'Oxygen pulls the shared electrons toward itself and leaves the hydrogens slightly positive. That imbalance is what makes water such a good solvent.',
  nh3: 'Nitrogen bonds three hydrogens and keeps one pair of electrons to itself, which pushes the molecule into a pyramid.',
  ch4: 'Carbon shares its four outer electrons with four hydrogens, which point to the corners of a tetrahedron.',
  hf: 'Fluorine pulls on electrons harder than any other element, so this bond is very lopsided.',
  n2: 'A triple bond, one of the strongest in chemistry. It is why the nitrogen that makes up most of the air hardly reacts.',
  co2: 'Two double bonds in a straight line. Each bond is lopsided, but they point in opposite directions and cancel out.',
  ch2o: 'A carbon–oxygen double bond, the group at the heart of sugars and many other molecules of life.',
  c2h4: 'A carbon–carbon double bond. Show the HOMO orbital to see the extra bond lying above and below the molecule.',
  c6h6: 'Six electrons shared around the whole ring instead of sitting in fixed bonds. It uses the quickest settings to stay smooth.',
  oh: 'A radical: one electron has no partner, which makes it extremely reactive.',
  o2: 'The oxygen we breathe has two unpaired electrons, which is why liquid oxygen sticks to a magnet.',
  glycine: 'The simplest amino acid, one of the twenty building blocks of every protein: an amine group and an acid group on one carbon.',
  pyridine: 'Benzene with one carbon swapped for a nitrogen. The nitrogen keeps a pair of electrons pointing out of the ring, which makes pyridine a base.',
  aspirin: 'Twenty-one atoms: a benzene ring carrying an acid group and an ester. Its highest filled orbitals are the ring’s π electrons.',
  glucose: 'The sugar our cells burn for energy, in the six-membered ring it forms in water. Twenty-four atoms, with five OH groups around the ring.',
  diamond:
    'Every carbon bonds to four others in a tetrahedron, without end. The ten carbons in the middle of one cubic cell are calculated; the bonds that would carry on into the crystal end in hydrogen atoms.',
  silicon:
    'The same structure as diamond with bigger atoms and longer, weaker bonds. Compare the gap between the highest filled and lowest empty levels with diamond’s: smaller, which is why silicon is a semiconductor.',
  graphite:
    'Carbon in flat sheets of hexagons, stacked loosely: pencil lead. The calculated flake has its π electrons spread over the whole sheet; show the HOMO to see them above and below the plane.',
  nacl: 'Sodium hands an electron to chlorine and the ions stack in a checkerboard of cubes. Eight ions are calculated; the cloud piles up on the chlorines, which hold the extra electrons.',
  mgo: 'The same checkerboard as table salt, with magnesium giving two electrons to each oxygen. The ions are smaller and more charged, so they pack tighter.',
};
