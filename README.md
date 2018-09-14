# Carbon nanotube and graphene sheet building interface

This package provide command line interface for building carbon nanotube
(TubeGen is required) and graphene sheet (CHARMM is required) with defects.

## Usage

Carbon nanotube

    ❯ python builder.py tubegen -h
    usage: builder.py tubegen [-h] [--tubegen TUBEGEN] [--chirality [N,N]]
                              [--prefix PREFIX] [--cap]

    optional arguments:
      -h, --help         show this help message and exit
      --tubegen TUBEGEN  path to tubegen executable
      --chirality [N,N]  chirality parameter n, m
      --prefix PREFIX    prefix for output file
      --cap              cap hydrogen if turned on

Graphene

    ❯ python builder.py graphene -h
    usage: builder.py graphene [-h] [--charmm CHARMM] [--nrings NRINGS] [--cap]
                               [--defect [0-1]] [--dense_defect]

    optional arguments:
      -h, --help       show this help message and exit
      --charmm CHARMM  path to CHARMM executable
      --nrings NRINGS  number of rings
      --cap            cap hydrogen if turned on
      --defect [0-1]   fraction of rings to have defect
      --dense_defect   favors defects to be near each other when turned on

## Example

    mkdir -p example/tube
    cd example/tube
    python ../../builder.py tubegen --tubegen ../../builder/tubegen-3.4/src/tubegen --chirality 5,2
    
![tube](example/tube/tube.png?raw=true "Tube")

    mkdir -p example/graphene
    cd example/graphene
    python ../../builder.py graphene --charmm ~/local/charmm/c40b1 --nrings 100 --defect 0
    
![graphene](example/graphene/graphene.png?raw=true "Graphene")


