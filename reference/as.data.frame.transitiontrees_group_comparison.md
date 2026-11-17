# Convert a group comparison to a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one of the three tables of a
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
result: the per-context tests, the omnibus tests, or the pairwise
distances between groups.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("pathways", "omnibus", "distance"),
  ...
)
```

## Arguments

- x:

  A `transitiontrees_group_comparison` object.

- row.names, optional:

  Ignored.

- what:

  `"pathways"` (default), `"omnibus"` or `"distance"`.

- ...:

  Ignored.

## Value

A data frame. For `what = "pathways"`, one row per tested context, with
the columns described in
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).
For `what = "omnibus"`, one row per axis with columns `axis`,
`statistic`, `value` and `p_value`. For `what = "distance"`, one row per
unordered pair of groups with columns `group_a`, `group_b` and
`distance` (the count-weighted symmetric Kullback-Leibler divergence of
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md),
nats).

## See also

[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
[`summary.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group_comparison.md).

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
cmp <- compare_groups(grp, iter = 20)
as.data.frame(cmp)
#>       pathway depth count_total count_High count_Low modal_High modal_Low flips
#> 1   synthesis     1         652        278       374  consensus consensus FALSE
#> 2     discuss     1        3951       2003      1948  consensus   discuss  TRUE
#> 3    cohesion     1        1695        938       757  consensus consensus FALSE
#> 4   consensus     1        6329       3392      2937       plan      plan FALSE
#> 5  coregulate     1        1970        894      1076       plan   discuss  TRUE
#> 6       adapt     1         509        141       368  consensus consensus FALSE
#> 7        plan     1        6157       2892      3265       plan      plan FALSE
#> 8     emotion     1        2837       1550      1287  consensus  cohesion  TRUE
#> 9     (start)     0       27533      13721     13812  consensus      plan  TRUE
#> 10    monitor     1        1433        633       800    discuss   discuss FALSE
#>       jsd_bits      jsd_p   jsd_padj  usage_g2    usage_p usage_padj
#> 1  0.102198705 0.04761905 0.05952381  13.86317 0.04761905 0.05357143
#> 2  0.065570921 0.04761905 0.05952381   1.42901 0.23809524 0.23809524
#> 3  0.042847830 0.04761905 0.05952381  22.14428 0.04761905 0.05357143
#> 4  0.016955946 0.04761905 0.05952381  47.94018 0.04761905 0.05357143
#> 5  0.016386725 0.04761905 0.05952381  16.86533 0.04761905 0.05357143
#> 6  0.015357190 0.09523810 0.10582011 105.30589 0.04761905 0.05357143
#> 7  0.013645832 0.04761905 0.05952381  26.39080 0.04761905 0.05357143
#> 8  0.013486858 0.04761905 0.05952381  29.61348 0.04761905 0.05357143
#> 9  0.007387220 0.04761905 0.05952381        NA         NA         NA
#> 10 0.001875704 0.95238095 0.95238095  19.42160 0.04761905 0.05357143
as.data.frame(cmp, what = "omnibus")
#>         axis                 statistic    value    p_value
#> 1 behavioral count-weighted JSD (bits) 874.1040 0.04761905
#> 2      usage                   sum G^2 282.9737 0.04761905
as.data.frame(cmp, what = "distance")
#>   group_a group_b   distance
#> 1    High     Low 0.04709345
```
