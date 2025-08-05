This is a companion repository for the preprint [arXiv:2504.18134
](https://arxiv.org/abs/2504.18134).

# Files and folders

`./CSPLS` contains the database of toric colorable seeds obtained in the article [The characterization of ($n−1$)-spheres with $n+4$ vertices having maximal Buchstaber number](https://doi.org/10.1515/crelle-2024-0027).

`./fan-giving maps` contains notebooks with the computations for the 58 (the boundary of a hexagon is ommited) seeds with Picard number $4$.

`./minimally non-fanlike` contains notebooks with the computations done for finding the minimally non-fanike seeds.

`./Characteristic_pair.py` and `./Simplicial_complex.py` are handmade Python libraries.

`./fanlike seeds` and `./minimally non-fanlike seeds` contain the list of facets of each fanlike and minimally non-fanlike seeds, respectively.

`./Oscar_fanlike_seeds`, respectively `./Oscar_minimally_non-fanlike_seeds` contains saved copies of the fanlike seeds, respectively minimally non-fanlike seeds, in Oscar format, theyr can be loaded un Julia with the command `load(<)`

`Oscar.ipynb` is a Julia notebook contaning all the complete non-singular fans based on a fanlike seed from the folder `./Oscar_fanlike_seeds`.
