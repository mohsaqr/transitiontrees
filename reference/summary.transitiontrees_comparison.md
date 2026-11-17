# Summarize a two-tree comparison

[`summary()`](https://rdrr.io/r/base/summary.html) reduces a
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
result to one row: the observed distance between the two trees, the
mean, standard deviation and 95th percentile of the permutation null
distribution, and the permutation p-value.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees_comparison` object.

- ...:

  Ignored.

## Value

An object of class `summary.transitiontrees_comparison`: a one-row data
frame with columns `distance`, `null_mean`, `null_sd`, `null_q95`,
`p_value`, `iter` and `n_pathways` (contexts in the per-context table).

## See also

[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
[`as.data.frame.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
summary(compare_trees(grp, iter = 20))
#> <transitiontrees_comparison summary>
#>  distance null_mean  null_sd null_q95 p_value iter n_pathways
#>    0.0471   0.00294 0.000556  0.00378  0.0476   20         10
```
