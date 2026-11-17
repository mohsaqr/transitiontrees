# Summarize the fit of a transition tree

`model_fit()` returns the fit statistics of a tree as one row:
log-likelihood, number of parameters, number of scored positions, AIC,
BIC and perplexity. The row puts trees with different depths, pruning or
smoothing on a common footing, and for a group of trees gives one row
per group.

Without `newdata` the statistics are in-sample. With `newdata` every
statistic is computed on those sequences, with AIC and BIC formed from
the held-out log-likelihood and the tree's number of parameters. The
number of parameters is the number of contexts times \\k - 1\\, for
\\k\\ states (see
[`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md)).

## Usage

``` r
model_fit(tree, newdata = NULL)
```

## Arguments

- tree:

  A `transitiontrees` or `transitiontrees_group` object.

- newdata:

  Sequences to evaluate the tree on. `NULL` (default) evaluates it
  in-sample.

## Value

A `data.frame` with columns `logLik`, `df`, `nobs`, `AIC`, `BIC` and
`perplexity`: one row for a tree, or one row per group, with a leading
`group` column, for a `transitiontrees_group`. When no position can be
scored, the statistics are `NA`.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or group,
or `newdata` has no usable sequence.

## See also

[`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md),
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
[`summary.transitiontrees_group()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group.md),
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
model_fit(tree)
#>      logLik df nobs      AIC      BIC perplexity
#> 1 -1511.774 76 1870 3175.548 3596.109   2.244394
model_fit(prune_tree(tree))
#>      logLik df nobs      AIC      BIC perplexity
#> 1 -1536.197 36 1870 3144.395 3343.607   2.273899
```
