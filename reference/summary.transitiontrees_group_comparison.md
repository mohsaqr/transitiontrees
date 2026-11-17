# Summarize a group comparison

[`summary()`](https://rdrr.io/r/base/summary.html) reduces a
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
result to one row: the omnibus behavioral and usage tests, and the
number of contexts whose adjusted p-value is below `alpha` on each axis
or whose modal next state differs between groups.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
summary(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `transitiontrees_group_comparison` object.

- alpha:

  False discovery rate cutoff used to count contexts. Default 0.05.

- ...:

  Ignored.

## Value

An object of class `summary.transitiontrees_group_comparison`: a one-row
data frame with columns `groups` (the group names joined by `" vs "`),
`n_contexts`, `iter`, `stratified`, `behavioral_jsd` and `behavioral_p`
(omnibus count-weighted Jensen-Shannon divergence, bits, and its
permutation p-value), `usage_g2` and `usage_p` (omnibus usage statistic
and its p-value), `alpha`, `n_behavioral` and `n_usage` (contexts with
adjusted p-value below `alpha`) and `n_flips` (contexts whose modal next
state differs). It prints as a short verdict.

## See also

[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
[`as.data.frame.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_group_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
summary(compare_groups(grp, iter = 20), alpha = 0.1)
#> Group comparison: High vs Low (10 contexts, 20 permutations)
#>   behavioral omnibus: JSD = 874.104 bits, p = 0.048
#>   usage omnibus     : G2 = 282.97, p = 0.048
#>   pathways with a behavioral difference (adjusted p < 0.10): 8
#>   pathways with a usage difference      (adjusted p < 0.10): 8
#>   pathways whose modal next state differs between groups: 4
```
