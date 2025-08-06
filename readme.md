This is a companion repository for the preprint [arXiv:2504.18134
](https://arxiv.org/abs/2504.18134).

# Files and folders

`./CSPLS` contains the database of toric colorable seeds obtained in the article [The characterization of ($n−1$)-spheres with $n+4$ vertices having maximal Buchstaber number](https://doi.org/10.1515/crelle-2024-0027).

`./fan-giving maps` contains Python notebooks with the computations for obtaining the 58 (the boundary of a hexagon is ommited) fanlike seeds with Picard number $4$.

`./minimally non-fanlike` contains Python notebooks with the computations performed for finding the minimally non-fanike seeds.

`Characteristic_pair.py` and `Simplicial_complex.py` are handmade Python libraries used in the two previous notebooks.

`fanlike seeds` and `minimally non-fanlike seeds` contain the list of facets of each fanlike seeds and minimally non-fanlike seeds, respectively.

`Oscar_fanlike_seeds` and `Oscar_minimally_non-fanlike_seeds` contain saved copies of the fanlike seeds and minimally non-fanlike seeds, respectively, in Oscar `SimplicialComplex`, one can regenerate these files by running the Julia script `generate_DB.jl`.

`fans_dim_<n>` are Julia notebooks containing the complete nonsingular fans obtained from the fanlike seeds. 

`neighbourly_smooth_polytope.ipynb` is a notebook containing the proof of existence of a smooth neighbourly polytope.

# Loading a seed

To load a seed from the saved files, use the following Julia code, after loading the Oscar package.
1. `K = load("./Oscar_fanlike_seeds/K_<i>_<n>")` for a fanlike seed where `i` is the number of the seed and `n` its dimension.
2. `L = load("./Oscar_minimally_non-fanlike_seeds/L_<i>_<n>")` for a fanlike seed where `i` is the number of the seed and `n` its dimension.


# Fans

To access the fans of dimension $n$, open the Julia notebook `fans_dim_<n>`.

For every fanlike seed with name `K_<i>_<n>`, its fans are encoded in a variable `fans_K_<i>_<n>`.

Fans whose characteristic map has a total of $k$ indeterminates are encoded as functions with $k$ input variables and output an Oscar `PolyhedralFan`.

### Example
`fans_K_6_3[1](3,8)` is the first fan of the fanlike seed `K_6_3` where the two indeterminates are set to $3$ and $8$.
