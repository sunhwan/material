import argparse
import sys

from .molecule import Molecule

def pdb2rtf(args):
    mol = Molecule.from_pdb(open(args.input))
    mol.to_rtf(open('output.rtf', 'w'))

def run(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    pconvert = subparsers.add_parser('pdb2rtf')
    pconvert.add_argument('input')
    pconvert.add_argument('--output', default="output")
    pconvert.set_defaults(func=pdb2rtf)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    run()
