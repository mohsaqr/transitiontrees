# Test where groups of sequences differ, context by context

`compare_groups()` tests where groups of sequences, such as high and low
achievers or two phases of a course, differ, context by context and for
the trees as a whole. Groups can differ in two distinct ways: they can
reach the same history and then continue differently (a behavioral
difference), or they can reach some histories more often than others do
(a usage difference). Both are tested for every context represented in
the groups' trees and for the trees as a whole, so that a difference can
be located in specific pathways rather than reported as one number.

The behavioral difference at a context is the Jensen-Shannon divergence
(Lin, 1991), in bits, between the groups' next-state distributions, each
group weighted by its count at that context. The usage difference is a
\\G^2\\ homogeneity statistic on the context's prevalence: its count
against the other windows of the same length in each group. The omnibus
statistics are the count-weighted sum of the divergences and the sum of
the \\G^2\\ values. All p-values come from one label permutation:
sequences are pooled, group labels are shuffled with group sizes fixed,
and every statistic is recomputed from raw counts. Per-context p-values
are adjusted for false discovery rate with the Benjamini-Hochberg (1995)
procedure, separately for each axis. When one person contributes several
sequences, a `block` identifier restricts the shuffle to within persons,
so the null distribution respects that nesting. The result is a
`transitiontrees_group_comparison` object whose tables are read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html);
[`summary()`](https://rdrr.io/r/base/summary.html) gives a one-row
overview and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
the largest behavioral differences.

`compare_groups()` handles two or more groups. For one overall distance
between two trees with a permutation test, use
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md);
it measures the symmetric Kullback-Leibler divergence, not the
Jensen-Shannon divergence used here.
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
maps the per-context differences of two groups next state by next state.

## Usage

``` r
compare_groups(group, iter = 999L, min_count = 1L, seed = 1L, block = NULL)
```

## Arguments

- group:

  A `transitiontrees_group` of unweighted trees, from
  `context_tree(..., group = )`.

- iter:

  Number of label permutations. Default 999.

- min_count:

  Contexts with a total count across groups below this are dropped from
  the table. Default 1.

- seed:

  Seed for the permutations, or `NULL` to use the current random-number
  state. Default 1. A non-`NULL` seed sets the session's random-number
  state.

- block:

  Block identifiers for a stratified permutation, one per sequence in
  group-then-row order. Default `NULL`: the block stored by
  `context_tree(..., block = )` is used when present, and labels are
  shuffled freely otherwise.

## Value

An object of class `transitiontrees_group_comparison`. Its tables are
returned by
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
`what = "pathways"` gives one row per tested context with columns
`pathway`, `depth`, `count_total`, one `count_<group>` and one
`modal_<group>` column per group (the most likely next state), `flips`
(whether the groups' modal next states differ), `jsd_bits`, `jsd_p`,
`jsd_padj`, `usage_g2`, `usage_p` and `usage_padj` (`NA` for the root,
which has no usage test), sorted by `jsd_bits`; `what = "omnibus"` gives
one row per axis; and `what = "distance"` gives one row per pair of
groups with their
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md).
[`summary()`](https://rdrr.io/r/base/summary.html) returns a one-row
overview and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
the largest behavioral differences.

## Conditions

- `transitiontrees_bad_input`:

  (error) `group` is not a `transitiontrees_group`, has fewer than two
  groups, or `block` does not have one entry per sequence.

- `transitiontrees_weighted_unsupported`:

  (error) a member tree was fitted with weights; the permutation
  recounts the raw sequences and would drop them.

- `transitiontrees_degenerate_permutation`:

  (warning) no block contains more than one group, so the stratified
  shuffle cannot move any label.

## References

Lin, J. (1991). Divergence measures based on the Shannon entropy. *IEEE
Transactions on Information Theory*, 37(1), 145-151.
[doi:10.1109/18.61115](https://doi.org/10.1109/18.61115)

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: A practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society: Series B (Methodological)*, 57(1),
289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## See also

[`summary.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group_comparison.md),
[`as.data.frame.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_group_comparison.md)
and
[`plot.transitiontrees_group_comparison`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_group_comparison.md)
for the result;
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for the pairwise test;
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
for the next-state difference map.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
cmp <- compare_groups(grp, iter = 20)
cmp
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
#>     pathway depth count_total modal_High modal_Low flips jsd_bits jsd_padj
#>   synthesis     1         652  consensus consensus FALSE    0.102    0.060
#>     discuss     1        3951  consensus   discuss  TRUE    0.066    0.060
#>    cohesion     1        1695  consensus consensus FALSE    0.043    0.060
#>   consensus     1        6329       plan      plan FALSE    0.017    0.060
#>  coregulate     1        1970       plan   discuss  TRUE    0.016    0.060
#>       adapt     1         509  consensus consensus FALSE    0.015    0.106
#>        plan     1        6157       plan      plan FALSE    0.014    0.060
#>     emotion     1        2837  consensus  cohesion  TRUE    0.013    0.060
#>     (start)     0       27533  consensus      plan  TRUE    0.007    0.060
#>     monitor     1        1433    discuss   discuss FALSE    0.002    0.952
#>  usage_padj
#>       0.054
#>       0.238
#>       0.054
#>       0.054
#>       0.054
#>       0.054
#>       0.054
#>       0.054
#>          NA
#>       0.054
summary(cmp)
#> Group comparison: High vs Low (10 contexts, 20 permutations)
#>   behavioral omnibus: JSD = 874.104 bits, p = 0.048
#>   usage omnibus     : G2 = 282.97, p = 0.048
#>   pathways with a behavioral difference (adjusted p < 0.05): 0
#>   pathways with a usage difference      (adjusted p < 0.05): 0
#>   pathways whose modal next state differs between groups: 4
```
