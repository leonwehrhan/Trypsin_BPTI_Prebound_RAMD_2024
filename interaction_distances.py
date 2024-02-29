import mdtraj as md
import numpy as np
import pandas as pd
import os
from utils import COMPLEXES, COMPLETED_TRJS


ATOM_IDX = {'TYR39-OH':300,
            'HIS40-O': 325,
            'PHE41-O': 345,
            'ASN97-O': 1154,
            'GLN194-NE2': 2512,
            'GLY195-N': 2517,
            'ASP196-N': 2524,
            'SER197-N': 2536,
            'SER197-OG': 2543,
            'SER212-O': 2745,
            'GLY214-N': 2770,
            'PRO13-O': 3422,
            'CYS14-O': 3432,
            'ETG15-N': 3433,
            'ETG15-O': 3445,
            'ARG17-N': 3456,
            'ARG17-NE': 3469,
            'ARG17-NH1': 3472,
            'ARG17-NH2': 3475,
            'ILE19-N':3499,
            'ARG39-NE': 3821,
            'ARG39-NH1': 3824,
            'ARG39-NH2': 3827}

PAIRS = [['SER212-O', 'ETG15-N'],
         ['SER197-OG', 'ETG15-N'],
         ['SER197-N', 'ETG15-O'],
         ['ASP196-N', 'ETG15-O'],
         ['GLY195-N', 'ETG15-O'],
         ['PHE41-O', 'ARG17-N'],
         ['GLY214-N', 'PRO13-O'],
         ['GLN194-NE2', 'CYS14-O'],
         ['TYR39-OH', 'ILE19-N']]

MIN_PAIRS = [[['ASN97-O', 'ARG39-NH1'],
              ['ASN97-O', 'ARG39-NH2'],
              ['ASN97-O', 'ARG39-NE']],
              [['HIS40-O', 'ARG17-NH1'],
               ['HIS40-O', 'ARG17-NH2'],
               ['HIS40-O', 'ARG17-NE']]]

Y151_RING_IDX = [1887, 1888, 1890, 1892, 1895, 1897]
R17_CZ_IDX = 3471

D_OUT='/home/lwehrhan/Documents/RAMD/interaction_distances'


def interaction_data(t):
    data = {}

    for p in PAIRS:
        s = f'{p[0]}--{p[1]}'
        d = dist(t, p[0], p[1])
        data[s] = d

    for p in MIN_PAIRS:
        s = f'{p[0][0]}--{p[0][1]}'
        d = min_dist(t, p)
        data[s] = d
    
    s = 'TYR151-s--ARG17-CZ'
    d = R17_Y151(t)
    data[s] = d
    
    return data


def R17_Y151(t):
    centroid = np.mean(t.xyz[:, Y151_RING_IDX], axis=1)
    cz = t.xyz[:, R17_CZ_IDX]

    d = np.linalg.norm(centroid - cz, axis=1)
    return d


def dist(t, atom_1, atom_2):

    idx_1 = ATOM_IDX[atom_1]
    idx_2 = ATOM_IDX[atom_2]

    d = md.compute_distances(t, [[idx_1, idx_2]])[:, 0]

    return d


def min_dist(t, pairs):

    idx_pairs = []
    for pair in pairs:
        idx_pair = [ATOM_IDX[x] for x in pair]
        idx_pairs.append(idx_pair)
    
    dists = md.compute_distances(t, idx_pairs)
    min_dists = np.min(dists, axis=1)

    return min_dists


if __name__ == '__main__':
    
    for X in COMPLEXES:
        for rep in COMPLETED_TRJS[X]:
            for i in COMPLETED_TRJS[X][rep]:
                t = md.load(os.path.join(f'TRYP_{X}_BPTI', f'{rep}_60', f'sim{i}', f'md{i}_aligned.xtc'), top=os.path.join(f'TRYP_{X}_BPTI', f'{rep}_60', 'traj.pdb'))
                di = interaction_data(t)
                df = pd.DataFrame(np.array([di[x] for x in di]).T, columns=[x for x in di])
                df.to_csv(os.path.join(D_OUT, f'{X}_{rep}_sim{i}.csv'), index=False)