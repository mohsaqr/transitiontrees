# Convert a tuning grid to a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the cross-validation grid, or only the selected configuration, as a
plain data frame.

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("grid", "best"),
  ...
)
```

## Arguments

- x:

  A `transitiontrees_tune` object.

- row.names, optional:

  Ignored.

- what:

  `"grid"` (default) for every configuration, or `"best"` for the
  selected one.

- ...:

  Ignored.

## Value

A data frame. For `what = "grid"`, one row per configuration, sorted by
held-out perplexity, with the columns described in
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).
For `what = "best"`, one row: the configuration with the lowest held-out
perplexity among those that scored on every fold, or zero rows when none
did.

## See also

[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`summary.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_tune.md).

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = 5,
                   prune = FALSE, folds = 3)
as.data.frame(tuned, what = "best")
#>   max_depth min_count                           smoothing prune    logLik
#> 1         2         5 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>   n_scored perplexity n_nodes_avg folds_failed
#> 1     1870   2.343845          13            0
```
