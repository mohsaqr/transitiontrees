# Count the observations a transition tree was fitted on

[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
observed positions in the sequences the tree was fitted on, the number
used by [`BIC()`](https://rdrr.io/r/stats/AIC.html). For a weighted fit
it is the weighted total. A held-out log-likelihood,
`logLik(tree, newdata)`, carries its own count of scored positions.

## Usage

``` r
# S3 method for class 'transitiontrees'
nobs(object, ...)
```

## Arguments

- object:

  A `transitiontrees` object.

- ...:

  Not used.

## Value

A single integer.

## See also

[`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md),
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
nobs(tree)
#> [1] 1870
```
