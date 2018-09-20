import argparse
import sys
import subprocess as sp
from math import cos, sin, radians
from functools import partial

sys.path.append('MoleculeX')

import numpy as np
from MoleculeX.moleculex import Molecule, Atom

_rtf_header = """
*  --------------------------------------------------------------------------  *
*          GRAPHENE                                                            *
*  --------------------------------------------------------------------------  *
*
36  1

MASS    23 HGR61    1.00800  ! aromatic H
MASS    61 CG2R61  12.01100  ! 6-mem aromatic C

DEFA FIRS NONE LAST NONE
AUTO ANGLES DIHE

"""

_prm_header = """*  --------------------------------------------------------------------------  *
*          GRAPHENE                                                            *
*  --------------------------------------------------------------------------  *
*

ATOMS
MASS    23 HGR61    1.00800  ! aromatic H
MASS    61 CG2R61  12.01100  ! 6-mem aromatic C
MASS    62 CG2R62  12.01100  ! 6-mem aromatic C

BONDS
CG2R61 CG2R61  305.00     1.3750 ! PROT benzene, JES 8/25/89
CG2R62 CG2R62  305.00     1.3750 ! PROT benzene, JES 8/25/89
CG2R61 CG2R62  305.00     1.3750 ! PROT benzene, JES 8/25/89
CG2R61 HGR61   340.00     1.0800 ! PROT phe,tyr JES 8/25/89

ANGLES
CG2R61 CG2R61 CG2R61   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R62 CG2R62 CG2R62   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R61 CG2R62 CG2R62   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R61 CG2R61 CG2R62   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R61 CG2R62 CG2R61   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R62 CG2R61 CG2R62   40.00    120.00   35.00   2.41620 ! PROT JES 8/25/89
CG2R61 CG2R61 HGR61    30.00    120.00   22.00   2.15250 ! PROT JES 8/25/89 benzene

DIHEDRALS
CG2R61 CG2R61 CG2R61 CG2R61     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R62 CG2R62 CG2R62 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R62 CG2R62 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R61 CG2R62 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R61 CG2R61 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R61 CG2R62 CG2R61     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R62 CG2R61 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R62 CG2R62 CG2R61     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R62 CG2R61 CG2R61 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R62 CG2R61 CG2R62 CG2R62     3.1000  2   180.00 ! PROT JES 8/25/89
CG2R61 CG2R61 CG2R61 HGR61      4.2000  2   180.00 ! PROT JES 8/25/89 benzene
HGR61  CG2R61 CG2R61 HGR61      2.4000  2   180.00 ! PROT JES 8/25/89 benzene

IMPROPER
CG2R61 X      X      HGR61     15.0000  0     0.00 ! PYRIDINE pyridine, yin

NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

HGR61    0.0       -0.0300     1.3582 ! benzene
CG2R61   0.0       -0.0700     1.9924 ! INDO/TRP
CG2R62   0.0       -0.0000     1.9924 ! INDO/TRP

END
"""

_tubegen_input = """
set format pdb
set chirality {chirality[0]},{chirality[1]}
set units angstrom
set gutter 1.6735,1.6735,0
set relax_tube yes
set cell_count 1,1,{cell_count_z}
generate
save {prefix}.pdb
exit
"""

_charmm_input = """* generate graphene sheet
*

ioformat ext
read rtf card name graphene.rtf
read para card name graphene.prm flex

read sequ gp 1
gene gp first none last none

read coor pdb name graphene.pdb resid

coor stat sele all end
set natom = ?nsel

hbuild sele type h* end

mini sd nstep 10000
mini abnr nstep 10000
write coor pdb name graphene.pdb

stop

energy
scalar type set 62 sele .not. ( bynu 1 .or. .bonded. bynu 1 ) end
define atom sele chem cg2r61 end

label loop

mini sd nstep 5000
mini abnr nstep 5000
mini sd nstep 5000
mini abnr nstep 5000
mini sd nstep 5000
mini abnr nstep 5000

scalar type set 61 sele ( ( .bonded. atom ) .and. .not. chem cg2r61 ) .subset. 1 end
define atom sele chem cg2r61 end

if ?nsel .lt. @natom goto loop

mini sd nstep 1000
mini abnr nstep 1000
mini sd nstep 1000
mini abnr nstep 1000

write coor pdb name graphene.pdb sele chem cg2r61 end

!write psf card name graphene.psf

stop
"""

def range_type(astr, min=0, max=101):
    value = float(astr)
    if min<= value <= max:
        return value
    else:
        raise argparse.ArgumentTypeError('value not in range %s-%s'%(min,max))

def csv_list_type(string):
    try:
        return [int(_) for _ in string.split(',')]
    except:
        raise ValueError('comma separated integers expected')

def to_rtf(mol, rtffile):
    v = []
    atom_names = {}
    count = {'C': 1, 'H': 1}
    chmtype_map = {'C': 'CG2R61', 'H': 'HGR61'}
    rtffile.write(_rtf_header)
    rtffile.write("RESI TUBE\n")
    for atom in mol.atoms():
        symbol = atom.element.symbol
        atom_names[atom.name] = ("%s%x" % (symbol, count[symbol])).upper()
        count[symbol] += 1
        chmtype = chmtype_map[symbol]
        rtffile.write("ATOM {} {} 0.0\n".format(atom_names[atom.name], chmtype))

    for atom1, atom2 in mol.bonds():
        rtffile.write("BOND {} {}\n".format(atom_names[atom1.name], atom_names[atom2.name]))
    rtffile.close()

def tubegen(args):
    pid = sp.Popen(args.tubegen, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, universal_newlines=True)
    (stdin, stdout) = (pid.stdin, pid.stdout)
    stdin.write(_tubegen_input.format(**args.__dict__))
    stdin.flush()
    print(stdout.read())

    pdbfile = '{}.pdb'.format(args.prefix)
    rtffile = '{}.rtf'.format(args.prefix)
    mol = Molecule.from_pdb(open(pdbfile))

    if args.cap:
        from math import sqrt
        atoms = list(mol.atoms())
        for atom in atoms:
            neighbors = [_ for _ in mol.graph.neighbors(atom) if _.element.symbol == 'C']
            if atom.element.symbol == 'C' and len(neighbors) == 2:
                r0 = np.array([atom.x, atom.y, atom.z])
                r1 = np.array([neighbors[0].x, neighbors[0].y, neighbors[0].z])
                r2 = np.array([neighbors[1].x, neighbors[1].y, neighbors[1].z])
                rc = (r1 + r2)/2
                vc = (r0 - rc) / np.linalg.norm(r0 - rc)
                rh = vc * 1.2 + r0

                cap = Atom(name='H', resname=atom.resname, resnr=atom.resnr, x=rh[0], y=rh[1], z=rh[2])
                mol.add_atom(cap)
                mol.add_bond(atom, cap)

            if atom.element.symbol == 'C' and len(neighbors) == 1:
                r0 = np.array([atom.x, atom.y, atom.z])
                neighbor = neighbors[0]
                r1 = np.array([neighbor.x, neighbor.y, neighbor.z])
                neighbor_neighbor = [_ for _ in mol.graph.neighbors(neighbor) if _.element.symbol == 'C' and _ != neighbor].pop()
                r2 = np.array([neighbor_neighbor.x, neighbor_neighbor.y, neighbor_neighbor.z])

                v1 = r0 - r1
                v2 = r1 - r2
                a = np.cross(v1, v2)
                c = cos(radians(60))
                s = sin(radians(60))
                ux, uy, uz = a
                cp = 1 - c
                R = np.array(((c + ux**2 * cp, ux * uy * cp - uz * s, ux * uz * cp + uy * s),
                              (uy * ux * cp + uz * s, c + uy**2 * cp, uy * uz * cp - ux * s),
                              (uz * ux * cp - uy * s, uz * uy * cp + ux * s, c + uz**2 * cp)))
                rh1 = np.dot(R, v1.T)
                rh1 = rh1 / np.linalg.norm(rh1) * 1.2

                a = v1
                ux, uy, uz = a
                c = cos(radians(120))
                s = sin(radians(120))
                cp = 1 - c
                R = np.array(((c + ux**2 * cp, ux * uy * cp - uz * s, ux * uz * cp + uy * s),
                              (uy * ux * cp + uz * s, c + uy**2 * cp, uy * uz * cp - ux * s),
                              (uz * ux * cp - uy * s, uz * uy * cp + ux * s, c + uz**2 * cp)))
                rh2 = np.dot(R, rh1.T)
                rh2 = rh2 / np.linalg.norm(rh2) * 1.2

                rh3 = np.dot(R, rh2.T)
                rh3 = rh3 / np.linalg.norm(rh3) * 1.2

                rh1 += r0
                rh2 += r0
                rh3 += r0

                cap = Atom(name='H', resname=atom.resname, resnr=atom.resnr, x=rh1[0], y=rh1[1], z=rh1[2])
                mol.add_atom(cap)
                mol.add_bond(atom, cap)

                cap = Atom(name='H', resname=atom.resname, resnr=atom.resnr, x=rh2[0], y=rh2[1], z=rh2[2])
                mol.add_atom(cap)
                mol.add_bond(atom, cap)

                cap = Atom(name='H', resname=atom.resname, resnr=atom.resnr, x=rh3[0], y=rh3[1], z=rh3[2])
                mol.add_atom(cap)
                mol.add_bond(atom, cap)

    to_rtf(mol, open(rtffile, 'w'))
    mol.to_pdb(open(pdbfile, 'w'))

def graphene(args):
    import graphene
    import networkx as nx
    from shutil import copy

    nrings = args.nrings
    defect = args.defect
    dense_defect = args.dense_defect
    cap_hydrogen = args.cap

    g = nx.Graph()
    g.defect_level = defect
    g.dense_defect = dense_defect
    graphene.add_unit(g)
    closed_node = []
    for i in range(nrings):
        node = graphene.find_neighbor(g)
        nodetype = g.node[node]['vertices']
        for j in range(nodetype):
            graphene.add_unit_neighbor(g, node)

        for n in g.nodes():
            if n not in closed_node and g.node[n]['vertices'] == g.degree(n):
                graphene.check_closure(g, n)
                closed_node.append(n)

    h = graphene.atom_graph(g)

    # cap hydrogen
    if cap_hydrogen:
        atoms = list(h.nodes())
        count = len(h.nodes())
        for n in atoms:
            if h.degree(n) != 3:
                nn = 3 - h.degree(n)
                for _ in range(nn):
                    h.add_node(count, type='H')
                    h.add_edge(n, count)
                    count += 1

    pos=nx.nx_pydot.graphviz_layout(h, prog='neato')
    pdbname = 'graphene.pdb'
    rtfname = 'graphene.rtf'
    atom_names = graphene.build_initial_pdb(g, pos, pdbname)
    graphene.build_topology(h, atom_names, rtfname)

    open('graphene.prm', 'w').write(_prm_header)
    open('mini.inp', 'w').write(_charmm_input)

    p = sp.Popen([args.charmm, '-i', 'mini.inp'])
    p.wait()

def run(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    ptubegen = subparsers.add_parser('tubegen')
    ptubegen.add_argument('--tubegen', default="./tubegen-3.4/src/tubegen", help="path to tubegen executable")
    ptubegen.add_argument('--chirality', default=[3, 3], type=csv_list_type, metavar='[N,N]', help="chirality parameter n, m")
    ptubegen.add_argument('--cell_count_z', default=1, type=int, metavar='N', help="cell count along Z")
    ptubegen.add_argument('--prefix', default="tube", help="prefix for output file")
    ptubegen.add_argument('--cap', default=False, action='store_true', help="cap hydrogen if turned on")
    ptubegen.set_defaults(func=tubegen)

    pgraphene = subparsers.add_parser('graphene')
    pgraphene.add_argument('--charmm', default="charmm", help="path to CHARMM executable")
    pgraphene.add_argument('--nrings', default=100, type=int, help="number of rings")
    pgraphene.add_argument('--cap', default=False, action='store_true', help="cap hydrogen if turned on")
    pgraphene.add_argument('--defect', default=0, type=partial(range_type, min=0, max=1), metavar='[0-1]', help="fraction of rings to have defect")
    pgraphene.add_argument('--dense_defect', default=False, action='store_true', help="favors defects to be near each other when turned on")
    pgraphene.set_defaults(func=graphene)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    run()
