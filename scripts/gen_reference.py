"""Generate basis-set JSON and PySCF reference fixtures.

Run with: python3 scripts/gen_reference.py   (needs `pip install pyscf`)

Basis functions are cartesian (6d), like the original Pople 6-31G* definition.
Every cartesian function is unit-normalized in this project, so matrices from
PySCF are rescaled by 1/sqrt(S_ii) before they are written out.
"""
import json
import os
import numpy as np
from pyscf import gto, scf, dft
from pyscf.dft import gen_grid

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
BASES = {'sto-3g': 'STO-3G', '6-31g': '6-31G', '6-31g*': '6-31G*'}


def write_basis():
    out = {}
    for key in BASES:
        out[key] = {}
        for el in ELEMENTS:
            shells = []
            for sh in gto.basis.load(key, el):
                l = sh[0]
                prims = np.array(sh[1:], dtype=float)
                exps = prims[:, 0]
                for col in range(1, prims.shape[1]):
                    coefs = prims[:, col]
                    keep = coefs != 0
                    shells.append({'l': l, 'exps': exps[keep].tolist(), 'coefs': coefs[keep].tolist()})
            out[key][el] = shells
    with open(os.path.join(ROOT, 'crates/qchem/data/basis-data.json'), 'w') as f:
        json.dump(out, f, separators=(',', ':'))


ANG = 1.0 / 0.52917721092
GEOMS = {
    'h2o': [('O', 0, 0, 0.1173), ('H', 0, 0.7572, -0.4692), ('H', 0, -0.7572, -0.4692)],
    'h2o_distorted': [('O', 0.05, 0, 0.11), ('H', 0, 0.80, -0.45), ('H', 0.1, -0.70, -0.50)],
    'ch4': [('C', 0, 0, 0), ('H', 0.6276, 0.6276, 0.6276), ('H', -0.6276, -0.6276, 0.6276),
            ('H', -0.6276, 0.6276, -0.6276), ('H', 0.6276, -0.6276, -0.70)],
    'nh3': [('N', 0, 0, 0.1162), ('H', 0, 0.9377, -0.2711), ('H', 0.8121, -0.4689, -0.2711),
            ('H', -0.8121, -0.4689, -0.2711)],
    'c2h4': [('C', 0, 0, 0.6695), ('C', 0, 0, -0.6695), ('H', 0, 0.9289, 1.2321), ('H', 0, -0.9289, 1.2321),
             ('H', 0, 0.9289, -1.2321), ('H', 0, -0.9289, -1.2321)],
    'n2': [('N', 0, 0, 0.56), ('N', 0, 0, -0.56)],
    'oh': [('O', 0, 0, 0), ('H', 0, 0, 0.98)],
    'ch2': [('C', 0, 0, 0.1), ('H', 0, 0.99, -0.45), ('H', 0, -0.99, -0.45)],
}


def make_mol(name, basis, spin=0):
    atoms = GEOMS[name]
    return gto.M(atom=[(a, (x, y, z)) for a, x, y, z in atoms], basis=basis, cart=True, spin=spin,
                 unit='Angstrom', verbose=0)


def geom_json(mol):
    return [{'symbol': mol.atom_symbol(i), 'xyz': mol.atom_coord(i).tolist()} for i in range(mol.natm)]


def normalized(mat, s):
    d = 1.0 / np.sqrt(np.diag(s))
    return (mat * d[:, None] * d[None, :])


def run_scf(mol, method, xc=None):
    if method == 'rhf':
        mf = scf.RHF(mol)
    elif method == 'uhf':
        mf = scf.UHF(mol)
    elif method == 'rks':
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.grids.level = 7
    mf.conv_tol = 1e-11
    mf.conv_tol_grad = 1e-7
    mf.kernel()
    assert mf.converged, (mol.atom, method, xc)
    return mf


def case(name, basis, method, xc=None, spin=0, grad=True, mats=False, eri=False):
    mol = make_mol(name, basis, spin)
    mf = run_scf(mol, method, xc)
    res = {
        'name': name, 'basis': basis, 'method': method, 'xc': xc, 'charge': 0, 'multiplicity': spin + 1,
        'atoms': geom_json(mol), 'nbf': mol.nao, 'energy': float(mf.e_tot),
        'dipole': mf.dip_moment(unit='AU', verbose=0).tolist(),
    }
    if method != 'uhf':
        res['moEnergies'] = mf.mo_energy.tolist()
    if grad:
        g = mf.nuc_grad_method()
        res['gradient'] = g.kernel().tolist()
    if mats:
        s = mol.intor('int1e_ovlp')
        res['S'] = normalized(s, s).tolist()
        res['T'] = normalized(mol.intor('int1e_kin'), s).tolist()
        res['V'] = normalized(mol.intor('int1e_nuc'), s).tolist()
    if eri:
        s = mol.intor('int1e_ovlp')
        d = 1.0 / np.sqrt(np.diag(s))
        e = mol.intor('int2e')
        e = np.einsum('ijkl,i,j,k,l->ijkl', e, d, d, d, d)
        flat = e.ravel()
        # Keep the fixture small: a fixed random sample of entries plus their flat indices.
        idx = np.sort(np.random.default_rng(7).choice(flat.size, size=min(flat.size, 4000), replace=False))
        res['eriIndex'] = idx.tolist()
        res['eri'] = flat[idx].tolist()
    return res


def write_fixtures():
    cases = [
        case('h2o', 'sto-3g', 'rhf', mats=True, eri=True),
        case('h2o', '6-31g*', 'rhf', mats=True, eri=True),
        case('h2o_distorted', '6-31g*', 'rhf'),
        case('ch4', '6-31g*', 'rhf'),
        case('nh3', '6-31g', 'rhf'),
        case('c2h4', 'sto-3g', 'rhf'),
        case('n2', '6-31g*', 'rhf'),
        case('oh', '6-31g', 'uhf', spin=1),
        case('ch2', '6-31g*', 'uhf', spin=2),
        case('h2o_distorted', '6-31g*', 'rks', xc='lda,vwn'),
        case('h2o_distorted', '6-31g*', 'rks', xc='pbe'),
        case('n2', '6-31g', 'rks', xc='pbe'),
        case('ch4', 'sto-3g', 'rks', xc='lda,vwn'),
    ]
    with open(os.path.join(ROOT, 'crates/qchem/tests/fixtures/pyscf-reference.json'), 'w') as f:
        json.dump(cases, f)


def write_lebedev():
    out = {}
    for n in (110, 194, 302, 590):
        g = gen_grid.MakeAngularGrid(n)
        out[str(n)] = [round(v, 15) for v in g.ravel().tolist()]  # x, y, z, w (weights sum to 1)
    with open(os.path.join(ROOT, 'crates/qchem/data/lebedev.json'), 'w') as f:
        json.dump(out, f, separators=(',', ':'))


if __name__ == '__main__':
    write_lebedev()
    write_basis()
    write_fixtures()
