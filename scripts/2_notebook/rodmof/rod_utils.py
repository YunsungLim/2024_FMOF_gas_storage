from pathlib import Path
import collections
from collections import defaultdict

import numpy as np

import ase
import ase.io
import ase.neighborlist

try:
    from ase.utils import natural_cutoffs
except Exception as e:
    from ase.neighborlist import natural_cutoffs

import pymatgen as mg
from pormake.framework import Framework


# Metal species.
METAL_LIKE = [
    "Li", "Be", "B", "Na", "Mg",
    "Al", "Si", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np",
    "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
]

def covalent_neigbor_list(
        atoms, scale=1.2, neglected_species=[], neglected_indices=[]):

    cutoffs = natural_cutoffs(atoms)
    cutoffs = [scale*c for c in cutoffs]
    # Remove radii to neglect them.
    species_indices = [
        i for i, a in enumerate(atoms) if a.symbol in neglected_species
    ]

    for i in neglected_indices+species_indices:
        cutoffs[i] = 0.0

    return ase.neighborlist.neighbor_list("ijD", atoms, cutoff=cutoffs)


def get_centroid(atoms):
    I, J, D = covalent_neigbor_list(atoms)

    nl = [[] for _ in atoms]
    for i, j, d in zip(I, J, D):
        nl[i].append((j, d))

    visited = [False for _ in atoms]
    q = collections.deque()

    abc_half = np.sum(atoms.get_cell(), axis=0) * 0.5

    positions = {}
    q.append((0, np.array([0.0, 0.0, 0.0])))
    while q:
        i, pos = q.pop()
        visited[i] = True
        positions[i] = pos
        for j, d in nl[i]:
            if not visited[j]:
                q.append((j, pos+d))
                visited[j] = True

    centroid = np.array([0.0, 0.0, 0.0])
    for v in positions.values():
        centroid += v
    centroid /= len(positions)
    first_position = atoms.get_positions()[0]
    return centroid + first_position


def clean_rod_mof(mof):
    MM_bonds = []
    # Find Metal-Metal Bonds
    MM_neighbors = defaultdict(list)
    new_bonds = []
    new_bonds_types = []
    
    symbols = mof.atoms.get_chemical_symbols()
    if any(['Lr' in symbols, 'No' in symbols]):
        for (a, b), t in zip(mof.bonds, mof.bond_types):
            if all([mof.atoms[a].symbol in ['No', 'Lr'], mof.atoms[b].symbol in ['No', 'Lr']]):
                MM_bonds.append([a,b])
            elif mof.atoms[a].symbol in ['No', 'Lr']:
                MM_neighbors[a].append([b,t])
            elif mof.atoms[b].symbol in ['No', 'Lr']:
                MM_neighbors[b].append([a,t])
            else:
                new_bonds.append((a,b))
                new_bonds_types.append(t)
    else:
        for (a, b), t in zip(mof.bonds, mof.bond_types):
            if all([mof.atoms[a].symbol in METAL_LIKE, mof.atoms[b].symbol in METAL_LIKE]):
                MM_bonds.append([a,b])
                
        for (a, b), t in zip(mof.bonds, mof.bond_types):
            if any([a in np.unique(MM_bonds), b in np.unique(MM_bonds)]):
                if all([mof.atoms[a].symbol in METAL_LIKE, mof.atoms[b].symbol in METAL_LIKE]):
                    continue
                elif mof.atoms[a].symbol in METAL_LIKE:
                    MM_neighbors[a].append([b,t])
                elif mof.atoms[b].symbol in METAL_LIKE:
                    MM_neighbors[b].append([a,t])
            else:
                new_bonds.append((a,b))
                new_bonds_types.append(t)


    for a, b in MM_bonds:
        mof.atoms[b].symbol = 'He'
        for i, t in MM_neighbors[b]:
            new_bonds.append((a, i))
            new_bonds_types.append(t)
        for i, t in MM_neighbors[a]:
            new_bonds.append((a, i))
            new_bonds_types.append(t)
            
    for a in mof.atoms:
        if a.symbol == 'No':
            a.symbol = 'Cl'
        elif a.symbol == 'Lr':
            a.symbol = 'O'

    new_mof = Framework(
        mof.atoms,
        new_bonds,
        new_bonds_types,
        mof.info,
        wrap=True
    )
    return new_mof
