# Print a summary of a tuning grid

[`print()`](https://rdrr.io/r/base/print.html) prints the number of
configurations and failed configurations, then the selected
configuration as a one-row table.

## Usage

``` r
# S3 method for class 'summary.transitiontrees_tune'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  A `summary.transitiontrees_tune` object.

- digits:

  Significant digits for numeric columns. Default 4.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_tune.md).

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:2, min_count = 5,
                   prune = FALSE, folds = 3)
print(summary(tuned), digits = 3)
#> <transitiontrees_tune summary>  2 configurations, 0 with a failed fold
#> selected (minimum held-out perplexity):
#>  max_depth min_count                           smoothing prune perplexity
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE       2.34
#>  n_nodes_avg
#>           13
```
