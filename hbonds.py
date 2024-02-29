import mdtraj as md
import numpy as np
import json
import os
from utils import mkdir
from utils import DPI, COMPLEXES, COLORS


D_OUT = ''


def hbonds_trj(t):
    hbonds = {}

    idx_tryp = t.top.select('resid 0 to 222')
    idx_bpti = t.top.select('resid 223 to 280')

    hb = md.wernet_nilsson(t)

    for frame in hb:
        for h in frame:
            if (h[0] in idx_tryp and h[2] in idx_bpti) or (h[0] in idx_bpti and h[2] in idx_tryp):
                a1 = t.top.atom(h[0])
                a2 = t.top.atom(h[2])

                if a1.is_sidechain:
                    s1 = f'{a1.residue.name}{a1.residue.resSeq}-s'
                else:
                    s1 = f'{a1.residue.name}{a1.residue.resSeq}-{a1.name}'
                
                if a2.is_sidechain:
                    s2 = f'{a2.residue.name}{a2.residue.resSeq}-s'
                else:
                    s2 = f'{a2.residue.name}{a2.residue.resSeq}-{a2.name}'
                
                s = f'{s1}--{s2}'

                if s in hbonds:
                    hbonds[s] += 1
                else:
                    hbonds[s] = 1
    
    hbonds = dict(sorted(hbonds.items(), key=lambda item: item[1], reverse=True))
    
    return hbonds


if __name__ == '__main__':

    for X in COMPLEXES:
        for i in range(20):

            t = md.load(os.path.join('fb', X, f'fully_bound_{i}', 'md_aligned.xtc'), top=os.path.join('fb', X, f'fully_bound_{i}', 'md_aligned.pdb'))
            hbonds = hbonds_trj(t)
            with open(os.path.join(D_OUT, f'fb_{X}_{i}.json'), 'w') as fp:
                json.dump(hbonds, fp)


            t = md.load(os.path.join('pb', X, f'pre_bound_{i}', 'md_aligned.xtc'), top=os.path.join('pb', X, f'pre_bound_{i}', 'md_aligned.pdb'))
            hbonds = hbonds_trj(t)
            with open(os.path.join(D_OUT, f'pb_{X}_{i}.json'), 'w') as fp:
                json.dump(hbonds, fp)
            
            print(X, i)


