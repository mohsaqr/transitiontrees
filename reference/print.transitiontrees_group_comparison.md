# Print a group comparison

[`print()`](https://rdrr.io/r/base/print.html) prints the groups and
permutation design, the omnibus behavioral and usage tests, and the `n`
contexts with the largest behavioral differences.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
print(x, n = 10L, digits = 3L, ...)
```

## Arguments

- x:

  A `transitiontrees_group_comparison` object.

- n:

  Number of contexts to print. Default 10.

- digits:

  Rounding digits for the printed tables. Default 3.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
[`as.data.frame.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_group_comparison.md)
for every context.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
print(compare_groups(grp, iter = 20), n = 3)
#> <transitiontrees_group_comparison>  2 groups, 20 permutations
#>   groups   : High, Low
#>   contexts : 10 tested
#> 
#> omnibus (permutation):
#>        axis                 statistic   value p_value
#>  behavioral count-weighted JSD (bits) 874.104   0.048
#>       usage                   sum G^2 282.974   0.048
#> 
#> top pathways by behavioral divergence (JSD):
#>    pathway depth count_total modal_High modal_Low flips jsd_bits jsd_padj
#>  synthesis     1         652  consensus consensus FALSE    0.102     0.06
#>    discuss     1        3951  consensus   discuss  TRUE    0.066     0.06
#>   cohesion     1        1695  consensus consensus FALSE    0.043     0.06
#>  usage_padj
#>       0.054
#>       0.238
#>       0.054
#> # ... 7 more (use as.data.frame(x) for the full table)
```
