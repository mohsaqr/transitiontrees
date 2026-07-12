# Number of Observations Used to Fit a context tree

Number of Observations Used to Fit a context tree

## Usage

``` r
# S3 method for class 'transitiontrees'
nobs(object, ...)
```

## Arguments

- object:

  A `transitiontrees`.

- ...:

  Ignored.

## Value

Integer. The number of state observations the tree was fitted on — the
(weight-adjusted) token total, equal to the `"nobs"` attribute of the
*in-sample*
[`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md).
(A held-out `logLik(newdata)` reports its own `"nobs"` — the positions
actually scored, e.g. excluding states outside the tree's alphabet —
which can be smaller.)
