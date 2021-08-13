import argparse
import sys
from pymatgen import Composition
from pymatgen.io.cif import CifParser, CifWriter
import re
import numpy as np
from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pymoo.model.problem import FunctionalProblem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.algorithms.genetic_algorithm import get_ages
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(
        description=(
            "Crystal Structure prediction by contact map,the coordination number of the cations and the ages of the individuals"
        )
    )

# data inputs
parser.add_argument(
    "--input",
    type=str,
    default="2-14-mp-236.input",
    metavar="PATH",
    help="input file",
)

parser.add_argument(
    "--template",
    default="",
    type=str,
    metavar="TEMPLATE",
    help="template cif file",
)

parser.add_argument(
    "--popsize",
    default=100,
    type=int,
    metavar="INT",
    help="popsize",
)

parser.add_argument(
    "--generation",
    default=1000,
    type=int,
    metavar="INT",
    help="generation",
)


parser.add_argument(
    "--random",
    default=50,
    type=int,
    metavar="INT",
    help="number of random structures",
)
args = parser.parse_args(sys.argv[1:])

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




def replace_sitexyz(text,coordinates):
    coordinates = coordinates.tolist()
    for i in range(len(coordinates)):
        coordinates[i]=str(coordinates[i])
    text1=text[text.rindex('_atom_site_occupancy')+len('_atom_site_occupancy')+3:]


    text1=text1.split('  ')
    text1=np.array(text1).reshape(int(len(text1)/7),7)
    h=0
    for i in range(len(text1)):
        if float(text1[i][3]) not in spe:
            text1[i][3] = coordinates[h]
            h+=1
        if float(text1[i][4]) not in spe:
            text1[i][4] = coordinates[h]
            h+=1
        if float(text1[i][5]) not in spe:
            text1[i][5] = coordinates[h]
            h+=1
    if h!=le[w]:
        print('error')
    text1=text1.flatten()
    text1='  '.join(text1)
    text=text.replace(text[text.rindex('_atom_site_occupancy')+len('_atom_site_occupancy')+3:],text1)
    return text


def cm_fitness(p):
    fitness=1111
    for z in p:
        if z in spe:
            fitness=9993
            break
    if fitness==1111:
        with open(files[w], "r") as fileobj:
            text = fileobj.read()
            text = replace_sitexyz(text, p)
            try:
                parser = CifParser.from_string(text)
                structure = parser.get_structures()[0]
                d2 = adj(structure)
            except ValueError:
                fitness=9990
            else:
                if d2.shape != ad[w].shape:
                    fitness = 9999
                else:
                    n = np.sum(d2 * ad[w] == 1)
                    n1 = np.sum(ad[w] == 1)
                    n2 = np.sum(d2 == 1)
                    try:
                        fitness = (-2) * n / (n1 + n2)
                    except ZeroDivisionError:
                        fitness = 9995
    return fitness



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

with open('spec/spec.txt','r') as fo:
    spe=fo.read()
    spe=re.findall(r'\d+.\d+', spe)
    spe = [float(x) for x in spe]

with open(args.input,'r') as fil:
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
latt = Lattice.from_para(lat[0], lat[1], lat[2], lat[3], lat[4], lat[5])

np.set_printoptions(threshold=np.inf)
fn=args.input.replace('.input','').replace('exp1/','')

if args.template!='':
    files=[args.template]
    print('template_file', args.template)
else:
    cmfit = []
    mps = []
    sh = []
    for k in range(args.random):
        C = pyxtal()
        C.from_random(3, spacegroup, elements, mols, lattice=latt)
        template_file = f"random_crystal/{fn}_template{k}.cif"
        random_crystal = f"random_crystal/{fn}_random{k}.cif"
        C.to_file(random_crystal)
        pmg_struc = C.to_pymatgen()
        W = CifWriter(pmg_struc, symprec=0.01)
        W.write_file(template_file)
        with open(template_file, 'r') as fl:
            tex = fl.read()
            tex = tex[tex.rindex('_atom_site_occupancy') + len('_atom_site_occupancy') + 3:]
            tex = tex.split('  ')
            tex = np.array(tex).reshape(int(len(tex) / 7), 7)
            sh.append(len(tex))
            el = []
            mp = []
            for s in range(len(tex)):
                el.append(tex[s][0])
                mp.append(float(tex[s][2]))
            wpel = []
            for e in elements:
                wpel.append(el.count(e))
            mpt = []
            for h in wpel:
                mpt.append(mp[:h])
                mp = mp[h:]
            mps.append(mpt)
            parser = CifParser(template_file)
            structure = parser.get_structures()[0]
            d1 = contactmap
            d2=adj(structure)
            if d1.shape != d2.shape:
                fitness = -9999
            else:
                n = np.sum(d2 * d1 == 1)
                n1 = np.sum(d1 == 1)
                n2 = np.sum(d2 == 1)

                try:
                    fitness = 2 * n / (n1 + n2)
                except ZeroDivisionError:
                    fitness = -9995
            cmfit.append(fitness)
    ind = []
    for z in range(len(mps)):
        tag = True
        for k in mps[z]:
            if k != sorted(k, reverse=True):
                tag = False
                break
        if tag == True:
            ind.append(z)

    sps = []
    cmfitne = []
    for a in ind:
        sps.append(sh[a])
        cmfitne.append(cmfit[a])

    index = []
    cmfitne2 = []
    for k in range(len(sps)):
        if sps[k] == min(sps):
            index.append(ind[k])
            cmfitne2.append(cmfitne[k])

    cmft=np.argsort(-np.array(cmfitne2))
    indtemp=[]
    for h in cmft:
        indtemp.append(index[h])

    index = indtemp[:5]
    print("finished generating and screening Random Crystal Structures")
    fbefore=[]
    files = []
    for z in index:
        template_file = f"random_crystal/{fn}_template{z}.cif"
        fbefore.append(cmfit[z])
        files.append(template_file)
    print('template_file', files)
    """
    print('contact fitness before', fbefore)
    bef=[]
    for s in files:
        bef.append(accu(args.input.replace('.input', '.cif'), s))
    print('contact fitness before', bef)
    """


ad=[]
le = []
co1 = []
co2 = []
s = 0
for i, file in enumerate(files):
    s = i
    ad.append(contactmap)
    with open(file, "r") as fileobj:
        text = fileobj.read()
        co = re.findall(r'\d+.\d+', text[text.rindex('occupancy'):])
        co1.append(co)
        cof = []
        cn=0
        for i in co:
            cof.append(float(i))
            if float(i) not in spe:
                cn+=1
        pf=cof[:]
        cof = np.array(cof).reshape(int(len(cof) / 3), 3)
        co2.append(cof)
        le.append(cn)


def coordination(p):
    fitness = 1111
    for z in p:
        if z in spe:
            fitness = 9993
            break
    if fitness == 1111:
        with open(files[w], "r") as fileobj:
            text = fileobj.read()
            text = replace_sitexyz(text, p)
            try:
                parser = CifParser.from_string(text)
                structure = parser.get_structures()[0]
                d2 = adj(structure)
            except ValueError:
                fitness = 9990
            else:
                if d2.shape != ad[w].shape:
                    fitness = 9999
                else:
                    with open(args.input.replace('.input', '.cc1'), 'r') as fo:
                        content = fo.read()
                        content=content[:len(content)-1]
                        fra = content.split('\n')
                        zsv=[]
                        for c in range(len(fra)):
                            if c>0:
                                v=fra[c].split()
                                zsv.append([v[0],v[2],v[3],v[4]])

                        fra=zsv[:]
                        ele = []
                        fracc = []
                        for s in fra:
                            ele.append(str(s[0]))
                            fracc.append([float(s[1]), float(s[2]), float(s[3])])

                        fra = np.array(fracc)
                        tap = CifParser(args.input.replace('.input', '.cif'))
                        tas = tap.get_structures()[0]
                        ain = []
                        aout = []
                        ine = []
                        oute = []
                        for m in range(len(fra)):
                            if 0<=fra[m][0]<round(lat[0], 5) and 0<=fra[m][1]<round(lat[1], 5) and 0<=fra[m][2]<round(lat[2], 5):
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
                    bl=int(sum(mols)/len(tas.cart_coords))
                    for t in range(len(d2)):
                        if np.sum(d2[t] == 1) > 0:
                            pw.append(np.sum(d2[t] == 1))
                            infra.append(structure.cart_coords[int(t/bl)])
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
                    pwnerr=0

                    for x in range(len(dist)):
                        cou = 0
                        for m in dist[x]:
                            if dyz[x][0] <= m <= dyz[x][1]:
                                cou += 1
                        if cou < pw[x]:
                            pwnerr+=pw[x]-cou
                        else:
                            pwnerr+=0


                    zz = '      '
                    if len(mols) == 3:
                        if elements[0] + zz[:6 - len(elements[1])] + elements[1] in yuzhi or elements[1] + zz[:6 - len(elements[0])] + elements[0] in yuzhi:
                            eqs = []
                            cz = [0, -1, 1]
                            for c in range(len(structure.sites)):
                                if structure.sites[c].as_dict()['species'][0]['element'] == elements[1]:
                                    for z in cz:
                                        for u in cz:
                                            for h in cz:
                                                if [z, u, h] != [0, 0, 0]:
                                                    eqs.append((structure.frac_coords[c] + [z, u, h]).tolist())

                            e1=elements[0]
                            z1=zz[:6 - len(elements[1])] + elements[1]
                            e2=elements[1]
                            z2=zz[:6 - len(elements[0])] + elements[0]
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
                    fitness=pwnerr
    return fitness



plo=[]
count1=0
def fitness_age(p):
    global count1
    count1=count1%args.popsize
    ages=get_ages()
    age_fitness=ages[count1]
    plo.append(age_fitness)
    count1+=1
    return age_fitness



best_co = []
best_fitness = []
c_fitness=[]
w = 0
outco=[]
pl=[]
predcif=[]
for i, file in enumerate(files):
    w = i
    problem = FunctionalProblem(le[i],
                                [cm_fitness, coordination,fitness_age],
                                xl=np.array([0.0] * le[i]),
                                xu=np.array([1.0] * le[i]),
                                parallelization='threads'
                                )

    algorithm = NSGA2(pop_size=args.popsize)

    resu = minimize(problem,
                   algorithm,
                   ('n_gen', args.generation),
                   verbose=False,
                    save_history=True)

    if isinstance(resu.X[0],float):
        resu.X=[resu.X]
    pl.append(len(resu.X))
    bes=[]
    for tx in resu.X:
        bes.append(cm_fitness(tx))

    bind=[]
    for z in range(len(bes)):
        if bes[z]==min(bes):
            bind.append(z)

    for v in range(len(bind)):
        best_x = resu.X[bind[v]]
        if cm_fitness(best_x) == 9993:
            c_a = cm_fitness(np.array(pf))
            c_fitness.append(c_a)
            best_co.append(pf)
            with open(files[w], "r") as fileobj:
                text = fileobj.read()
            pr=f'_predicted{v}.cif'
            res = fileobj.name.replace('random_crystal', 'resultcif').replace('.cif', pr)
            predcif.append(res)
            fp = open(res, "w")
            fp.truncate()
            with open(res, "w") as fw:
                fw.write(text)
        else:
            c_a = cm_fitness(best_x)
            c_fitness.append(c_a)

            with open(files[w], "r") as fileobj:
                text = fileobj.read()
                text = replace_sitexyz(text, best_x)
            best_xval = re.findall(r'\d+.\d+', text[text.rindex('occupancy'):])
            best_xval = [float(e) for e in best_xval]
            best_co.append(best_xval)
            pr = f'_predicted{v}.cif'
            res = fileobj.name.replace('random_crystal', 'resultcif').replace('.cif', pr)
            predcif.append(res)
            fp = open(res, "w")
            fp.truncate()
            with open(res, "w") as fw:
                fw.write(text)

for n in predcif:
    print(n)


