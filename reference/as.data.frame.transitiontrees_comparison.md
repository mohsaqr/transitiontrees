# Convert a two-tree comparison to a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the per-context divergence table or the permutation null distribution of
a
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
result.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("pathways", "null"),
  ...
)
```

## Arguments

- x:

  A `transitiontrees_comparison` object.

- row.names, optional:

  Ignored.

- what:

  `"pathways"` (default) for the per-context table, or `"null"` for the
  permuted distances.

- ...:

  Ignored.

## Value

A data frame. For `what = "pathways"`, one row per context present in
either tree, sorted by `divergence_sym` (largest first), with columns
`pathway`, `count_a` and `count_b` (the context's count in each tree, 0
when absent), `divergence_ab` and `divergence_ba` (the two directed
Kullback-Leibler divergences, nats) and `divergence_sym` (their mean).
For `what = "null"`, one row per permutation with columns `permutation`
and `distance`.

## See also

[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
[`summary.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
cmp <- compare_trees(grp, iter = 20)
as.data.frame(cmp)
#>       pathway count_a count_b divergence_ab divergence_ba divergence_sym
#> 1   synthesis     278     374   0.300094275   0.310806923    0.305450599
#> 2     discuss    2003    1948   0.177769139   0.204033457    0.190901298
#> 3    cohesion     938     757   0.119512707   0.114971189    0.117241948
#> 4       adapt     141     368   0.043562551   0.094339054    0.068950802
#> 5   consensus    3392    2937   0.049933054   0.045829733    0.047881394
#> 6  coregulate     894    1076   0.045521559   0.047136936    0.046329248
#> 7        plan    2892    3265   0.041779504   0.035907622    0.038843563
#> 8     emotion    1550    1287   0.034618429   0.034624558    0.034621494
#> 9     (start)   13721   13812   0.019555937   0.022035405    0.020795671
#> 10    monitor     633     800   0.005243522   0.005331551    0.005287536
head(as.data.frame(cmp, what = "null"))
#>   permutation    distance
#> 1           1 0.002055196
#> 2           2 0.003808779
#> 3           3 0.003212274
#> 4           4 0.003508735
#> 5           5 0.002851955
#> 6           6 0.002360542
```
