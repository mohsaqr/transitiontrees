# Print a summary of a two-tree comparison

[`print()`](https://rdrr.io/r/base/print.html) prints the one-row
summary of a
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
result as a table.

## Usage

``` r
# S3 method for class 'summary.transitiontrees_comparison'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `summary.transitiontrees_comparison` object.

- digits:

  Significant digits for numeric columns. Default 3.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
print(summary(compare_trees(grp, iter = 20)), digits = 2)
#> <transitiontrees_comparison summary>
#>  distance null_mean null_sd null_q95 p_value iter n_pathways
#>     0.047    0.0029 0.00056   0.0038   0.048   20         10
```
