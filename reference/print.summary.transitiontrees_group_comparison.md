# Print a summary of a group comparison

[`print()`](https://rdrr.io/r/base/print.html) prints the one-row
summary of a
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
result as a short verdict: the omnibus tests and the counts of contexts
that differ.

## Usage

``` r
# S3 method for class 'summary.transitiontrees_group_comparison'
print(x, ...)
```

## Arguments

- x:

  A `summary.transitiontrees_group_comparison` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
print(summary(compare_groups(grp, iter = 20)))
#> Group comparison: High vs Low (10 contexts, 20 permutations)
#>   behavioral omnibus: JSD = 874.104 bits, p = 0.048
#>   usage omnibus     : G2 = 282.97, p = 0.048
#>   pathways with a behavioral difference (adjusted p < 0.05): 0
#>   pathways with a usage difference      (adjusted p < 0.05): 0
#>   pathways whose modal next state differs between groups: 4
```
