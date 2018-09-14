import pytest
import moleculex as mx

def test_atom():
    v = mx.Atom(symbol='H')
    assert v.element.number == 1
    assert v.element.symbol == 'H'
    assert v.element.mass == 1.00800

    v = mx.Atom(number=1)
    assert v.element.number == 1
    assert v.element.symbol == 'H'
    assert v.element.mass == 1.00800

def test_atom_element_number_mismatch():
    with pytest.raises(AssertionError):
        v = mx.Atom(number=1, symbol='C')

def test_atom_unknown_atom():
    with pytest.raises(KeyError):
        v = mx.Atom(symbol='X')