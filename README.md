# Galois-Groups-of-cubic-surfaces-supporting-data
This repository contains the supporting data and code for the paper [TO ADD].


## Code

### Monodromy computations

This part of the code is implemented in `SageMath` and relies of the [`lefschetz-family` package](https://github.com/ericpipha/lefschetz-family), in which general methods not specific to cubic and quartic surfaces have been implemented.in the `Fibration` class.

There is an example of a symmetric cubic computation, code for computing the monodromy of all $G$-invariant cubic surfaces for $G$ a subgroup of $S_4$.

There is an example of a crystallographic quartic computation, where we have also included the period matrix of the base fibre to shorten the computation.
To recompute the period matrix (which takes around one hour), delete the file `fibre_32.9.sobj`.

### Discriminant computations

This part of the code is implemented in `julia` and relies on OSCAR *TODO: add link*.

*TODO explain what is invluded and how to use it*

## Data

### Symmetric cubic surfaces

*TODO: describe data*

### Crystallographic cubic surfaces

*TODO: describe data*

### Symmetric quartic surfaces

Instead of including the symmetric quartic surfaces, we refer to [Table 7 of the paper]() *TODO add link to page 29* to match to the corresponding crystallographic group.

### Crystallographic quartic surfaces

*TODO: describe data*