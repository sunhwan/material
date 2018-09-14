elements_number = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
}

elements_symbol = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O'
}

elements_mass = {
    1: 1.00800,
    6: 12.01100,
    7: 14.00700,
    8: 15.99900,
}

class Element:
    """Element"""

    symbol = None
    number = None
    mass = None

    def __init__(self, **kwargs):
        if 'symbol' in kwargs:
            self.symbol = kwargs['symbol']
        if 'number' in kwargs:
            self.number = kwargs['number']

        assert self.symbol or self.number

        if not self.symbol:
            self.symbol = elements_symbol[self.number]
        if not self.number:
            self.number = elements_number[self.symbol]

        assert elements_symbol[self.number] == self.symbol
        assert elements_number[self.symbol] == self.number

        self.mass = elements_mass[self.number]

    def __eq__(self, other):
        if isinstance(other, Element):
            return self.symbol == other.symbol
        raise NotImplementedError
