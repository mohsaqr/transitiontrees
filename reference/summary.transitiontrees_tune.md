# Summarize a tuning grid

[`summary()`](https://rdrr.io/r/base/summary.html) reduces a
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
grid to one row: the number of configurations scored, the number with a
failed cross-validation fold, and the settings, held-out perplexity and
mean size of the selected configuration.

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees_tune` object.

- ...:

  Ignored.

## Value

An object of class `summary.transitiontrees_tune`: a one-row data frame
with columns `n_configurations`, `n_failed`, `max_depth`, `min_count`,
`smoothing`, `prune`, `perplexity` and `n_nodes_avg`. The last six
describe the selected configuration and are `NA` when no configuration
scored on every fold.

## See also

[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`as.data.frame.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_tune.md).

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = 5,
                   prune = c(FALSE, TRUE), folds = 3)
summary(tuned)
#> <transitiontrees_tune summary>  6 configurations, 0 with a failed fold
#> selected (minimum held-out perplexity):
#>  max_depth min_count                           smoothing prune perplexity
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE      2.344
#>  n_nodes_avg
#>           13
```
