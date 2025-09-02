# Galois-Groups-of-cubic-surfaces-supporting-data

This repository contains the supporting data and code for the paper *TODO: add link*.


## Code

Code that can be used to reproduce the results of the paper is included in the `code` repository, and we give here a brief overview of how it may be used.

### Monodromy computations

This part of the code is implemented in `SageMath` and relies on the [`lefschetz-family` package](https://github.com/ericpipha/lefschetz-family), in which general methods not specific to cubic and quartic surfaces have been implemented.in the `Fibration` class.

The file `cubic_example.sage` contains an example of a symmetric cubic computation, followed by code for computing the monodromy of all $G$-invariant cubic surfaces for $G$ a subgroup of $S_4$.

There is an example of a crystallographic quartic computation, where we have also included the period matrix of the base fibre to shorten the computation.
To recompute the period matrix (which takes around one hour), delete the file `fibre_32.9.sobj`.
The same piece of code can be used to compute the monodromy of all crystallographic families of quartic surfaces, by running, e.g., 
```
sage quartic_example.sage "19.2"
```
from inside the `sage_code` repository (passing no argument runs `32.9` by default).
All cases except `2.1`, `3.1`, `8.1`, `10.1`, `26.1` and `27.1` only involve integration of operators of order 9 or lower, and should go through in a few hours.
The remaining cases require lengthier computation.

### Discriminant computations

This part of the code is implemented in `julia` and relies on [OSCAR](https://www.oscar-system.org). The package UnirationalDiscriminants.jl implements methods for computing the discriminant of a linear system of projective hypersurfaces. 

To use the package, make sure you are in the folder UnirationalDiscriminants.
You can check this by accessing the shell in a julia terminal by clicking ";" and then using `pwd`:

```shell
pwd
```

If this does not show the path to UnirationalDiscriminants, move to the correct folder by using `cd`:
```shell
cd path/to/folder/UnirationalDiscriminants
```

Hit backspace to exit the shell. Next, use the following commands:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using UnirationalDiscriminants
```

Some functionalities of the package are illustrated in the file `examples.jl`.


## Data

### Symmetric cubic surfaces

This folder contains computational results for the families of symmetric cubic surfaces appearing in [Table 1 of the paper]() *TODO add link to page 3*. There is one txt file for each group in that table, and such a file contains
- a list of generators of the group as a subgroup of $\mathrm{GL}_4(\mathbb{Z})$,
- the basis we chose for the linear system of invariants,
- the permutation of the 27 lines induced by each of the generators of the group, 
- the irreducible components of the discriminant hypersurface,
- generators of the Galois group as a subgroup of $S_{27}$. 
In addition, the file `discriminants.txt` contains the output of the discriminant computation. The file `multidegs.txt` contains the multigradings and multidegrees used for the computation of the component ${\rm pr}_2(Y_0)$ of the discriminant, in cases where this is a hypersurface. 

### Crystallographic cubic surfaces

This folder contains computational results for the families of crystallographic cubic surfaces appearing in [Table 3 of the paper]() *TODO add link to page 29*. There is one txt file for each group in that table, and such a file contains
- a list of generators of the group as a subgroup of $\mathrm{GL}_4(\mathbb{Z})$,
- the basis we chose for the linear system of invariants,
- a list of equivalent crystallographic groups,
- the permutation of the 27 lines induced by each of the generators of the group, 
- the irreducible components of the discriminant hypersurface,
- generators of the Galois group as a subgroup of $S_{27}$. 
In addition, the file `discriminants.txt` contains the output of the discriminant computation when the crystallographic group is not equivalent to a symmetric group. In the latter case, the discriminant is that for the corresponding symmetric group up to a linear change of coordinates.  

### Symmetric quartic surfaces

Instead of including the symmetric quartic surfaces, we refer to [Table 7 of the paper]() *TODO add link to page 29* to match to the corresponding crystallographic group. 

### Crystallographic quartic surfaces

This folder contains computational results for the families of crystallographic quartic surfaces appearing in [Table 8 of the paper]() *TODO add link to page 29*. There is one txt file for each group in that table, and such a file contains
- a list of generators of the group as a subgroup of $\mathrm{GL}_4(\mathbb{Z})$,
- the basis we chose for the linear system of invariants,
- a list of equivalent crystallographic groups,
- the intersection product of the homology of the base fibre,
- a basis for the generic Picard lattice,
- the homology class of the hyperplane section,
- the lattice isomorphism induced by each of the generators of the group, 
- a list of homology classes of all conics contained in the generic Picard lattice, 
- the irreducible components of the discriminant hypersurface 
- if we were not able to compute the discriminant, we include the restriction of the discriminant to a randomly chosen line, 
- generators of the monodromy group on the second homology group $H_2(X)$,
- the induced generators on the monodromy on the generic Picard lattice,
- a list of generators of the monodromy group on the Picard lattice, represented as permutations on the listed conics. 

In addition, the file `discriminants.txt` contains the output of the discriminant computation, in cases where we managed to compute it. The file `discriminants_lines.txt` contains details of the discriminant computations after restricting to a randomly chosen line in the linear system. 

