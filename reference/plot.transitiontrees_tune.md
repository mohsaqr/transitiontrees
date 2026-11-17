# Plot the cross-validated perplexity of a tuning grid

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
held-out perplexity surface of a
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
grid: perplexity against maximum depth, one line per minimum count, with
one panel per smoothing method (and per pruning choice when both were
tried). The selected configuration is marked with a star. The plot shows
where added depth stops lowering held-out perplexity and how sensitive
that point is to the other settings.

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
plot(x, ...)
```

## Arguments

- x:

  A `transitiontrees_tune` object from
  [`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

- ...:

  Ignored.

## Value

A ggplot object: held-out perplexity (y) against `max_depth` (x),
colored by `min_count`, faceted by smoothing and pruning.

## See also

[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`summary.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_tune.md).

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = c(5, 10),
                   prune = FALSE, folds = 3)
plot(tuned)
```
