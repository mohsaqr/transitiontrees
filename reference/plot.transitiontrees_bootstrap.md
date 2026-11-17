# Plot bootstrap intervals of pathway informativeness

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a forest
plot of each pathway's likelihood-ratio statistic \\G^2\\ against its
parent context: the bootstrap mean with its percentile interval, one row
per pathway, colored and ordered by whether the pathway is stable,
informative, both or neither. A dashed line marks the chi-square
critical value used by the bootstrap; a pathway whose interval lies
wholly to its right departs reproducibly from its shorter history.

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
plot(x, top = 25L, min_stability = NULL, ...)
```

## Arguments

- x:

  A `transitiontrees_bootstrap` object.

- top:

  Maximum number of pathways shown. Default 25.

- min_stability:

  Smallest `stability_rate` shown. Default `NULL`: the
  `stability_threshold` of the bootstrap.

- ...:

  Ignored.

## Value

A ggplot object: one row per pathway with the mean \\G^2\\ and its
interval.

## Conditions

- `transitiontrees_empty_result`:

  (error) no pathway with a \\G^2\\ reaches `min_stability`.

## See also

[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
[`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md)
for the resample distributions.

## Examples

``` r
boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
                           iter = 50)
plot(boot, min_stability = 0.5)
```
