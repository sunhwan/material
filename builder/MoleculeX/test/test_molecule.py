import pytest
import moleculex as mx

def test_add_atom():
    m = mx.Molecule()
    v = mx.Atom()
    m.add_atom(v)
    assert len(m.atoms()) == 1

def test_add_bond():
    m = mx.Molecule()
    v = mx.Atom()
    w = mx.Atom()
    m.add_atom(v)
    m.add_atom(w)
    m.add_bond(v, w)
    assert len(m.bonds()) == 1
    assert len(m.atoms()) == 2
