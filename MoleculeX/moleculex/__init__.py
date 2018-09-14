"""
MoleculeX

>>> import moleculex as mx
>>> G = mx.Molecule()
>>> G.add_bond(1,2)
>>> G.add_atom(42)
>>> print(sorted(G.atoms()))
[1, 2, 42]
>>> print(sorted(G.bonds()))
[(1, 2)]
"""

from .molecule import Molecule
from .atom import Atom
from .element import Element
