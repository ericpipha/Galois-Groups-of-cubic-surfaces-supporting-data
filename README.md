# Galois-Groups-of-cubic-surfaces-supporting-data

This repository contains the supporting data and code for the paper *TODO: add link*.


## Code

Code that can be used to reproduce the results of the paper is included in the `code` repository, and we give here a brief overview of how it may be used.

### Monodromy computations

This part of the code is implemented in `SageMath` and relies on the [`lefschetz-family` package](https://github.com/ericpipha/lefschetz-family), in which general methods not specific to cubic and quartic surfaces have been implemented.in the `Fibration` class.

There is an example of a symmetric cubic computation, code for computing the monodromy of all $G$-invariant cubic surfaces for $G$ a subgroup of $S_4$.

There is an example of a crystallographic quartic computation, where we have also included the period matrix of the base fibre to shorten the computation.
To recompute the period matrix (which takes around one hour), delete the file `fibre_32.9.sobj`.
The same piece of code can be used to compute the monodromy of all crystallographic families of quartic surfaces, by changing the value of `key` at the top of the file.
All cases except `2.1`, `3.1`, `8.1`, `10.1`, `26.1` and `27.1` only involve integration of operators of order 9 or lower, and should go through in a few hours.
The remaining cases require lengthier computation.

### Discriminant computations

This part of the code is implemented in `julia` and relies on OSCAR *TODO: add link*.

*TODO explain what is included and how to use it*

## Data

### Symmetric cubic surfaces

*TODO: describe data*

### Crystallographic cubic surfaces

*TODO: describe data*

### Symmetric quartic surfaces

Instead of including the symmetric quartic surfaces, we refer to [Table 7 of the paper]() *TODO add link to page 29* to match to the corresponding crystallographic group.

### Crystallographic quartic surfaces

*TODO: describe data*


The cases `2.1` and `3.1` will be added shortly, the computation of the full monodromy representation on $H_2(X)$ being still in progress.