# Plot the permutation distribution of a two-tree comparison

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
permutation null distribution of the distance between two context trees
as a histogram, with the observed distance as a vertical line and the
permutation p-value in the subtitle. An observed line far in the right
tail indicates that the two groups predict differently beyond what
random reassignment of sequences produces.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
plot(x, bins = 30L, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison` object from
  [`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md).

- bins:

  Number of histogram bins. Default 30.

- ...:

  Ignored.

## Value

A ggplot object: a histogram of the permuted distances with the observed
distance marked.

## See also

[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
[`as.data.frame.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_comparison.md)
for the permuted values.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
plot(compare_trees(grp, iter = 20), bins = 10)
```
