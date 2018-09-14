import networkx as nx

from .atom import Atom

class Molecule():
    """Molecule class"""

    def __init__(self):
        self.graph = nx.Graph()

    def add_atom(self, n):
        """Add a single atom"""
        self.graph.add_node(n)

    def add_bond(self, v, w):
        """Add a single bond"""
        self.graph.add_edge(v, w)

    def atoms(self):
        """Return a list of atoms"""
        return self.graph.nodes()

    def bonds(self):
        """Return a list of bonds"""
        return self.graph.edges()

    @staticmethod
    def from_pdb(pdbfile):
        molecule = None
        in_chain = False
        current_chain = None
        current_residue = None
        for line in pdbfile:
            token = line[:6].strip()
            if not in_chain and (token in ('ATOM', 'HETATM')):
                in_chain = True
                if molecule is None:
                    molecule = Molecule()
                    atomDict = {}
                else:
                    raise NotImplementedError('Parsing multiple molecule in a single PDB file is not implemented yet.')
            if in_chain and token == 'TER':
                in_chain = False

            if in_chain and (token in ('ATOM', 'HETATM')):
                number = int(line[6:11])
                name = line[11:16].strip()
                altloc = line[16:17].strip()
                resname = line[17:20].strip()
                chainid = line[21:22]
                resnr = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occupancy = float(line[54:60])
                bfactor = float(line[60:66])
                segid = line[72:76].strip()

                if altloc and altloc != 'A':
                    continue

                atom = Atom(name=name, resname=resname, resnr=resnr, x=x, y=y, z=z, occupancy=occupancy, bfactor=bfactor)
                atomDict[number] = atom
                molecule.add_atom(atom)

            if token == 'CONECT':
                entries = [int(_) for _ in line.strip().split()[1:]]
                atom1 = atomDict[entries[0]]
                for i in range(1, len(entries)):
                    atom2 = atomDict[entries[i]]
                    molecule.add_bond(atom1, atom2)

        return molecule

    def to_rtf(self, rtffile):
        raise NotImplementedError

    def to_pdb(self, pdbfile):
        atomnr = 0
        _atom = "ATOM  {atomnr:5}{name:^6}{resname:<4.4}{chainid}{resnr:4}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}\n"
        for atom in self.atoms():
            resname = atom.resname
            atomnr += 1
            chainid = 'A'
            name = atom.name
            resnr = atom.resnr
            x = atom.x
            y = atom.y
            z = atom.z
            try:
                occupancy = atom.occupancy
                bfactor = atom.bfactor
            except:
                occupancy = 0
                bfactor = 0
            pdbfile.write(_atom.format(**locals()))
            #ATOM      1  N   MET A   1     109.412 102.640  60.093  1.00 73.05           N

