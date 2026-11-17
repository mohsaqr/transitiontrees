# Print a tuning grid

[`print()`](https://rdrr.io/r/base/print.html) prints the number of
configurations, the first `n` rows of the grid (lowest held-out
perplexity first) and the selected configuration.

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `transitiontrees_tune` object.

- n:

  Number of configurations to print. Default 10.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`as.data.frame.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_tune.md)
for the full grid.

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:2, min_count = 5,
                   prune = FALSE, folds = 3)
print(tuned, n = 3)
#> <transitiontrees_tune>  2 configurations
#>  max_depth min_count                           smoothing prune    logLik
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>          1         5 floor(ymin=0.001, rule=interpolate) FALSE -1659.904
#>  n_scored perplexity n_nodes_avg folds_failed
#>      1870   2.343845          13            0
#>      1870   2.429412           4            0
#> 
#> best (min perplexity):
#>  max_depth min_count                           smoothing prune    logLik
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>  n_scored perplexity n_nodes_avg folds_failed
#>      1870   2.343845          13            0
```
