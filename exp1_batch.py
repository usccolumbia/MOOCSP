import os
from natsort import natsorted
import numpy as np

with open('style.ini', "r") as fileobj:
    yuzhi = fileobj.read()

def all_distances(coords1, coords2):
    """
    Returns the distances between two lists of coordinates

    Args:
        coords1: First set of cartesian coordinates.
        coords2: Second set of cartesian coordinates.

    Returns:
        2d array of cartesian distances. E.g the distance between
        coords1[i] and coords2[j] is distances[i,j]
    """
    c1 = np.array(coords1)
    c2 = np.array(coords2)
    z = (c1[:, None, :] - c2[None, :, :]) ** 2
    return np.sum(z, axis=-1) ** 0.5



def adj(structure):
    d2=all_distances(structure.cart_coords,structure.cart_coords)
    nShortBonds = 0
    for i, s1 in enumerate(structure.sites):
        for j, s2 in enumerate(structure.sites):
            if i == j:
                d2[i, j] = 0
            else:
                e1 = s1.as_dict()['species'][0]['element']
                e2 = s2.as_dict()['species'][0]['element']
                z='      '
                z1=z[:6-len(e2)]+e2
                z2=z[:6-len(e1)]+e1
                if e1 + z1 in yuzhi:
                    max_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 15:yuzhi.rindex(
                        e1 + z1) + len(e1 + z1) + 22])
                    min_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 4:yuzhi.rindex(
                        e1 + z1) + len(e1 + z1) + 11])
                elif e2 + z2 in yuzhi:
                    max_len = float(yuzhi[
                                    yuzhi.rindex(e2 + z2) + len(e2 + z2) + 15:yuzhi.rindex(
                                        e2 + z2) + len(e2 + z2) + 22])
                    min_len = float(yuzhi[
                                    yuzhi.rindex(e2 + z2) + len(e2 + z2) + 4:yuzhi.rindex(
                                        e2 + z2) + len(e2 + z2) + 11])
                else:
                    max_len = 0.00000
                    min_len = 0.00000
                if min_len <= d2[i, j] <= max_len:
                    d2[i, j] = 1
                elif d2[i,j]>max_len:
                    d2[i, j] = 0
                else:
                    nShortBonds+=1
    if nShortBonds>=1:
        d2 = all_distances(structure.cart_coords,structure.cart_coords).flatten()
    d2 = d2.astype(int)
    return d2


for root,dirs,filess in os.walk("exp1"):
    filess = natsorted(filess)
files=[]
for i in filess:
    if '.cif' in i:
        te = root + '/' + i
        if i.replace('.cif','.input') not in filess:
            os.system(f'python cm_test.py --cif {te}')
        k1 = te.replace('.cif', '.input')
        os.system(f'python age_fitness.py --input {k1}')
        for root2, dirs2, filess2 in os.walk("resultcif"):
            filess2 = natsorted(filess2)
        for z in filess2:
            if i.replace('.cif','') in z and 'measure' not in z:
                z1=root2+'/'+z
                os.system(f'python measure.py --cif {te} --predicted {z1}')
        files.append(te)
