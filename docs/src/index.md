# MovingWeightedLeastSquares.jl

## Overview

MovingWeightedLeastSquares.jl is a package that provides an implementation of the moving weighted least squares method.
A (very) nice short description of this method by Andy Nealen can be found [here](http://nealen.de/projects/mls/asapmls.pdf).

Let $\theta(d): \mathbb{R}^+ \rightarrow \mathbb{R}^+$ be a weighting function of the method.
Very often the there will be an $\varepsilon \in \mathbb{R}$, such that $\forall \delta > \varepsilon: \theta(\delta) = 0$.
Whenever we use $\varepsilon$ or `EPS` in this document, we mean the cutoff distance for the weighting function.
An example of a good weighting function is $\theta(d) = \exp(d^2 / a^2)$, where $a$ is the average distance between sample input data.

Interaction with this package is done mostly via structure `MwlsObject` and its subclasses `MwlsKdObject`, `MwlsCllObject` and `MwlsNaiveObject`.
This interface is similar to the interface of "interpolations objects" from [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).

The difference between the subclasses is the solution of the range search problem.
`MwlsKdObject` solves the range search problem by using a k-d tree created by Kristoffer Carlsson, see [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl).
`MwlsCllObject` solves the range search problem by using a cell linked list, which is implemented in this package.
If the cell linked list is needed 
`MwlsNaiveObject` solves the range search problem naively.  
**TL;DR**: use anything but `MwlsNaiveObject`.

## Installation

At the moment this packaged can be installed by manually cloning the package

```
Pkg.clone("https://github.com/vutunganh/MovingWeightedLeastSquares.jl")
```
