import argparse
import sys
import subprocess as sp
from functools import partial

sys.path.append('MoleculeX')

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
set cell_count {cell_count[0]},{cell_count[1]},{cell_count[2]}
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
                r0 = [atom.x, atom.y, atom.z]
                r1 = [neighbors[0].x, neighbors[0].y, neighbors[0].z]
                r2 = [neighbors[1].x, neighbors[1].y, neighbors[1].z]
                rc = [(r1[i] + r2[i])/2 for i in range(3)] # (r1 + r2)/2
                vc = [(r0[i] - rc[i]) for i in range(3)] # r0 - rc
                norm = sqrt(vc[0]**2 + vc[1]**2 + vc[2]**2)
                vc = [vc[i] / norm for i in range(3)] # np.linalg.norm(vc)
                rh = [(vc[i] * 1.2 + r0[i]) for i in range(3)] # vc * 1.2 + r0

                cap = Atom(name='H', resname=atom.resname, resnr=atom.resnr, x=rh[0], y=rh[1], z=rh[2])
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
    ptubegen.add_argument('--tubegen', default="./tubegen-3.4/src/tubegen")
    ptubegen.add_argument('--chirality', default=[3, 3], type=csv_list_type)
    ptubegen.add_argument('--cell_count', default=[1, 1, 1], type=csv_list_type)
    ptubegen.add_argument('--prefix', default="tube")
    ptubegen.add_argument('--cap', default=False, action='store_true')
    ptubegen.set_defaults(func=tubegen)

    pgraphene = subparsers.add_parser('graphene')
    pgraphene.add_argument('--charmm', default="charmm")
    pgraphene.add_argument('--nrings', default=100, type=int)
    pgraphene.add_argument('--cap', default=False, action='store_true')
    pgraphene.add_argument('--defect', default=0, type=partial(range_type, min=0, max=1), metavar='[0-1]')
    pgraphene.add_argument('--dense_defect', default=False, action='store_true')
    pgraphene.set_defaults(func=graphene)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    run()
