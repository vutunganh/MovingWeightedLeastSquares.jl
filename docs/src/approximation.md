# Approximation

Let `obj` be a subclass of `MwlsObject`.
The function call `obj(pt::Point)` returns the approximated function value at `pt`.
This calls the [`approximate`](@ref) function.

## Example

Let's initialize a dataset with input data `xs` and output data `fs`.

```@example approx
using MovingWeightedLeastSquares # hide
xs = collect(-2:0.1:2);
fs = [sin(x) for x in xs];
```

Now let's construct an approximation object `obj`.
Let's choose a weighting function ``$\theta(d) = \exp(d^2)$``.

```@example approx
obj = mwlsKd(xs, fs, 0.5, (d, e) -> (exp(-d^2)));
```

The approximation at 1 can be obtained by using

```@example approx
obj(1)
```

If a different range of neighbor data is needed, then use

```@example
obj(1; dist = 1)
```

## Relevant documentation

```@docs
calcMwlsCoefficients
approximate
```

# Approximation of derivative

Let `obj` be a subclass of `MwlsObject`.
To calculate an approximation at `pt` `obj(pt::Point)` should be used.

