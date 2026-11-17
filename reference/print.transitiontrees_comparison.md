# Print a two-tree comparison

[`print()`](https://rdrr.io/r/base/print.html) prints the number of
permutations, the observed distance, the mean of the null distribution,
the p-value, and the `n` contexts that contribute most to the distance.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
print(x, digits = 3L, n = 6L, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison` object.

- digits:

  Significant digits for the distances and rounding digits for the
  per-context table. Default 3.

- n:

  Number of contexts to print. Default 6.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
[`as.data.frame.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_comparison.md)
for every context.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
print(compare_trees(grp, iter = 20), n = 3)
#> <transitiontrees_comparison>  iter = 20
#>   observed distance : 0.0471
#>   null mean         : 0.00294
#>   p-value           : 0.0476
#> 
#> top divergent pathways:
#>    pathway count_a count_b divergence_ab divergence_ba divergence_sym
#>  synthesis     278     374         0.300         0.311          0.305
#>    discuss    2003    1948         0.178         0.204          0.191
#>   cohesion     938     757         0.120         0.115          0.117
```
