import sys
import argparse
import numpy as np
import re
from pymatgen import Composition
from pymatgen.io.cif import CifParser
from cifparser import CifParserExpand
from sklearn.metrics import mean_squared_error, mean_absolute_error
from pymatgen.analysis.structure_matcher import StructureMatcher as SM


parser = argparse.ArgumentParser(
        description=(
            "measure predicted structure"
            "calculate contact map accuracy,coordination error,rmsd,mae and rms"
        )
    )

parser.add_argument(
    "--cif",
    type=str,
    default="2-14-mp-236.cif",
    metavar="PATH",
    help="Path to target cif file",
)

parser.add_argument(
    "--predicted",
    type=str,
    default="2-14-mp-236_predicted.cif",
    metavar="PATH",
    help="Path to predicted cif file",
)

args = parser.parse_args(sys.argv[1:])

"""
# Usage:
# python measure.py --cif 2-14-mp-236.cif --predicted 2-14-mp-236_predicted.cif

"""


with open('style.ini', "r") as fileobj:
    yuzhi = fileobj.read()


with open('spec/spec.txt','r') as fo:
    spe=fo.read()
    spe=re.findall(r'\d+.\d+', spe)
    spe = [float(x) for x in spe]


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


def coordination(pred):
    with open(pred, "r") as fileobj:
        text = fileobj.read()
        try:
            parser = CifParser.from_string(text)
            structure = parser.get_structures()[0]
            d2 = adj(structure)
        except ValueError:
            fitness = 9990
        else:
            if d2.shape != contactmap.shape:
                fitness = 9999
            else:
                with open(args.cif.replace('.cif', '.cc1'), 'r') as fo:
                    content = fo.read()
                    content = content[:len(content) - 1]
                    fra = content.split('\n')
                    zsv = []
                    for c in range(len(fra)):
                        if c > 0:
                            v = fra[c].split()
                            zsv.append([v[0], v[2], v[3], v[4]])

                    fra = zsv[:]
                    ele = []
                    fracc = []
                    for s in fra:
                        ele.append(str(s[0]))
                        fracc.append([float(s[1]), float(s[2]), float(s[3])])

                    fra = np.array(fracc)
                    tap = CifParser(args.cif)
                    tas = tap.get_structures()[0]
                    ain = []
                    aout = []
                    ine = []
                    oute = []
                    for m in range(len(fra)):
                        if 0 <= fra[m][0] < round(lat[0], 5) and 0 <= fra[m][1] < round(lat[1], 5) and 0 <= fra[m][
                            2] < round(lat[2], 5):
                            ain.append(fra[m])
                            ine.append(ele[m])
                        else:
                            aout.append(fra[m])
                            oute.append(ele[m])

                d2 = all_distances(ain, aout)
                nShortBonds = 0
                for i, s1 in enumerate(ain):
                    for j, s2 in enumerate(aout):
                        e1 = ine[i]
                        e2 = oute[j]
                        z = '      '
                        z1 = z[:6 - len(e2)] + e2
                        z2 = z[:6 - len(e1)] + e1
                        if e1 + z1 in yuzhi:
                            max_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 15:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 22])
                            min_len_v = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 4:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 11])
                        elif e2 + z2 in yuzhi:
                            max_len = float(yuzhi[
                                            yuzhi.rindex(e2 + z2) + len(e2 + z2) + 15:yuzhi.rindex(
                                                e2 + z2) + len(e2 + z2) + 22])
                            min_len_v = float(yuzhi[
                                              yuzhi.rindex(e2 + z2) + len(e2 + z2) + 4:yuzhi.rindex(
                                                  e2 + z2) + len(e2 + z2) + 11])
                        else:
                            max_len = 0.00000
                            min_len_v = 0.00000

                        min_len = min_len_v
                        if min_len <= d2[i, j] <= max_len:
                            d2[i, j] = 1
                        elif d2[i, j] > max_len:
                            d2[i, j] = 0
                        else:
                            nShortBonds += 1

                pw = []
                infra = []
                dyz = []
                bl = int(sum(mols) / len(tas.cart_coords))
                for t in range(len(d2)):
                    if np.sum(d2[t] == 1) > 0:
                        pw.append(np.sum(d2[t] == 1))
                        infra.append(structure.cart_coords[int(t / bl)])
                        e1 = structure.sites[t].as_dict()['species'][0]['element']
                        e2 = elements[-1]
                        z = '      '
                        z1 = z[:6 - len(e2)] + e2
                        z2 = z[:6 - len(e1)] + e1
                        if e1 + z1 in yuzhi:
                            max_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 15:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 22])
                            min_len_v = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 4:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 11])
                        elif e2 + z2 in yuzhi:
                            max_len = float(yuzhi[
                                            yuzhi.rindex(e2 + z2) + len(e2 + z2) + 15:yuzhi.rindex(
                                                e2 + z2) + len(e2 + z2) + 22])
                            min_len_v = float(yuzhi[
                                              yuzhi.rindex(e2 + z2) + len(e2 + z2) + 4:yuzhi.rindex(
                                                  e2 + z2) + len(e2 + z2) + 11])
                        else:
                            max_len = 0.00000
                            min_len_v = 0.00000
                        dyz.append([min_len_v, max_len])

                eqs = []
                cz = [0, -1, 1]
                for c in range(len(structure.sites)):
                    if structure.sites[c].as_dict()['species'][0]['element'] == elements[-1]:
                        for z in cz:
                            for u in cz:
                                for h in cz:
                                    if [z, u, h] != [0, 0, 0]:
                                        eqs.append((structure.frac_coords[c] + [z, u, h]).tolist())

                eqs_bk = []
                [eqs_bk.append(i) for i in eqs if i not in eqs_bk]
                eqs = eqs_bk[:]
                eqscar = structure.lattice.get_cartesian_coords(eqs)
                dist = all_distances(infra, eqscar)
                pwnerr = 0

                for x in range(len(dist)):
                    cou = 0
                    for m in dist[x]:
                        if dyz[x][0] <= m <= dyz[x][1]:
                            cou += 1
                    if cou < pw[x]:
                        pwnerr += pw[x] - cou
                    else:
                        pwnerr += 0

                zz = '      '
                if len(mols) == 3:
                    if elements[0] + zz[:6 - len(elements[1])] + elements[1] in yuzhi or elements[1] + zz[:6 - len(
                            elements[0])] + elements[0] in yuzhi:
                        eqs = []
                        cz = [0, -1, 1]
                        for c in range(len(structure.sites)):
                            if structure.sites[c].as_dict()['species'][0]['element'] == elements[1]:
                                for z in cz:
                                    for u in cz:
                                        for h in cz:
                                            if [z, u, h] != [0, 0, 0]:
                                                eqs.append((structure.frac_coords[c] + [z, u, h]).tolist())

                        e1 = elements[0]
                        z1 = zz[:6 - len(elements[1])] + elements[1]
                        e2 = elements[1]
                        z2 = zz[:6 - len(elements[0])] + elements[0]
                        if e1 + z1 in yuzhi:
                            max_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 15:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 22])
                            min_len_v = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 4:yuzhi.rindex(
                                e1 + z1) + len(e1 + z1) + 11])
                        elif e2 + z2 in yuzhi:
                            max_len = float(yuzhi[
                                            yuzhi.rindex(e2 + z2) + len(e2 + z2) + 15:yuzhi.rindex(
                                                e2 + z2) + len(e2 + z2) + 22])
                            min_len_v = float(yuzhi[
                                              yuzhi.rindex(e2 + z2) + len(e2 + z2) + 4:yuzhi.rindex(
                                                  e2 + z2) + len(e2 + z2) + 11])

                        eqs_bk = []
                        [eqs_bk.append(i) for i in eqs if i not in eqs_bk]
                        eqs = eqs_bk[:]
                        eqscar = structure.lattice.get_cartesian_coords(eqs)
                        dist = all_distances(infra[:int(mols[0] / bl)], eqscar)

                        for x in range(len(dist)):
                            for m in dist[x]:
                                if min_len_v <= m <= max_len:
                                    pwnerr = 9997
                                    break
                fitness = pwnerr
    return fitness


def accu(cif,other):
    with open(cif, "r") as fileobj:
        text1 = fileobj.read()
        co = re.findall(r'\d+\.\d+', text1[text1.rindex('occupancy'):])
        cof1 = []
        for i in co:
            cof1.append(float(i))
        cof1 = np.array(cof1).reshape(int(len(cof1) / 3), 3)
        zs = text1[text1.rindex('_atom_site_occupancy') + len('_atom_site_occupancy') + 3:]
        mpl = zs.split('  ')
        mpl = np.array(mpl).reshape(int(len(mpl) / 7), 7)
        mulpl = []
        for u in range(len(mpl)):
            mulpl.append(mpl[u][2])

    with open(other, "r") as fileobj:
        textpr = fileobj.read()
        copr = re.findall(r'\d+\.\d+', textpr[textpr.rindex('occupancy'):])
        coprf1 = []
        for i in copr:
            coprf1.append(float(i))
        pwerr=coordination(other)
        coprf1 = np.array(coprf1).reshape(int(len(coprf1) / 3), 3)


    parser1 = CifParserExpand(cif)
    structure1 = parser1.get_structures()[0]

    try:
        parser2 = CifParserExpand(other)
        structure2 = parser2.get_structures()[0]
        d2=adj(structure2)
    except ValueError:
        fitness = -9990
    else:
        d1 = contactmap
        if d2.shape != d1.shape:
            fitness = -9999
        else:
            n = np.sum(d2 * d1 == 1)
            n1 = np.sum(d1 == 1)
            n2 = np.sum(d2 == 1)

            try:
                fitness = 2 * n / (n1 + n2)
            except ZeroDivisionError:
                fitness = -9995


    sm = SM(ltol=0.6, stol=0.6, angle_tol=20)
    try:
        rms = sm.get_rms_dist(structure1, structure2)[0]
    except TypeError:
        rms = 'None'

    if len(cof1) == len(coprf1):
        cos2 = parser2.cos
        cos2t = []
        for i in range(len(cos2)):
            cos2t.append(cos2[-i - 1])
        tem = []
        for k in cos2t:
            tem += k

        te = np.array(tem).flatten()
        fz = []
        for w in mulpl:
            fz.append(te[:3 * int(w)])
            te = te[3 * int(w):]
        te = fz
        cos2f = []
        for i in range(len(te)):
            cos2f.append(te[i].reshape(int(len(te[i]) / 3), 3))

        co3 = []
        for i in range(len(cof1)):
            for j in range(len(cos2f[i])):
                co3.append(np.sqrt(mean_squared_error(cof1[i], cos2f[i][j])))

        tempe = []
        for w in mulpl:
            tempe.append(co3[:int(w)])
            co3 = co3[int(w):]
        co3 = tempe

        co4 = []
        for u in range(len(co3)):
            k = np.argmin(co3[u])
            co4.append(cos2f[u][k])

        min_rmse = np.sqrt(mean_squared_error(cof1.flatten(), np.array(co4).flatten()))
        min_mae = mean_absolute_error(cof1.flatten(), np.array(co4).flatten())
    else:
        min_rmse=5555
        min_mae=5555
    return fitness,pwerr,min_rmse,min_mae,rms

with open(args.cif.replace('.cif','.input'),'r') as fil:
    inp=fil.read()
    formula=inp[inp.rindex('formula:')+len('formula:'):inp.rindex('spacegroup')].replace('\n','')
    spacegroup=int(re.findall('\d+',inp[inp.rindex('spacegroup'):inp.rindex('contactmap')])[0])
    lattice = re.findall(r'\d+\.\d+', inp[inp.rindex('lattice-abc'):inp.rindex('contactmap:')])
    dim=re.findall('\d+',inp[inp.rindex('('):inp.rindex(')')])
    contactmap=re.findall('\d+',inp[inp.rindex(')'):])
    for i in range(len(contactmap)):
        contactmap[i]=int(contactmap[i])
    contactmap=np.array(contactmap).reshape(int(dim[0]),int(dim[0]))



comp=Composition(formula)
elements=list(comp.as_dict().keys())
mols= list(comp.as_dict().values())
mols=[int(x) for x in mols]


lat=[]
for k in lattice:
    lat.append(float(k))

fitness,coordination_error,min_rmse,min_mae,rms=accu(args.cif,args.predicted)

fp = open(args.predicted.replace('.cif','_measure'), "w")
fp.truncate()
with open(args.predicted.replace('.cif','_measure'), "w") as fw:
    fw.write(args.predicted+'\n')
    fw.write('contact map accuracy:' + str(fitness) + '\n')
    fw.write('coordination error:' + str(coordination_error) + '\n')
    fw.write('rmsd:'+str(min_rmse)+'\n')
    fw.write('mae:'+str(min_mae)+'\n')
    fw.write('rms:' + str(rms) + '\n')


