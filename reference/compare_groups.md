# Compare Groups of Sequences for Structural Differences

Test how a set of fitted group trees (a `transitiontrees_group` from
`context_tree(..., group =)`) differ, on two complementary axes:

- behavioral:

  given the *same* context, do the groups predict a different next
  state? Measured per context by the count-weighted Jensen-Shannon
  divergence (bits) across the groups' next-state distributions.

- usage:

  do the groups *reach* a context at different rates? Measured per
  context by a \\G^2\\ homogeneity statistic on its prevalence (its
  share of each group's positions).

Significance is assessed by a label-permutation null throughout
(per-pathway and omnibus), with Benjamini-Hochberg FDR on the
per-pathway p-values.

## Usage

``` r
compare_groups(group, iter = 999L, min_count = 1L, seed = 1L, block = NULL)
```

## Arguments

- group:

  A `transitiontrees_group` of unweighted trees (a weighted member is
  refused: the permutation recounts the raw sequences).

- iter:

  Integer. Number of label permutations. Default 999.

- min_count:

  Integer. Drop contexts whose total count across all groups is below
  this. Default 1.

- seed:

  Integer or `NULL`. RNG seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html); `NULL` uses the
  current RNG state. Default 1.

- block:

  Block ids for a **stratified** permutation: group labels are shuffled
  only *within* each block, so the null respects nested /
  repeated-measures structure (e.g. several sequences from one subject)
  and holds any between-block difference fixed. Normally you do not pass
  this — fit with `context_tree(..., block = )` and it is carried on the
  object and used automatically. Passing a vector here (one id per
  sequence, in pooled group-then-row order) overrides that. `NULL` with
  no stored block shuffles labels freely.

## Value

A `transitiontrees_group_comparison`: a list with

- pathways:

  Per-context data.frame sorted by `jsd_bits` descending, with columns
  `pathway`, `depth`, `count_total`, one `count_<group>` and one
  `modal_<group>` column per group (most likely next state, ties broken
  by alphabet order), `flips` (do the groups' modal next states
  disagree?), `jsd_bits`, `jsd_p`, `jsd_padj`, `usage_g2`, `usage_p`,
  `usage_padj`. `usage_*` is `NA` for the root, which has no prevalence
  test.

- omnibus:

  Two-row data.frame (one row per axis) with columns `axis`,
  `statistic`, `value` and `p_value` (permutation p-value: one plus the
  number of permuted statistics at least the observed one, over
  `iter + 1`).

- distance_matrix:

  K x K symmetric-KL distance matrix between the groups (from
  [`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)).

- groups, iter, seed, n_contexts, stratified:

  Configuration; `stratified` is `TRUE` when a `block` was used.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `pathways` table; the print method shows the omnibus table and the
top pathways; [`summary()`](https://rdrr.io/r/base/summary.html) prints
a compact verdict.

## Details

The permutation pools every sequence, shuffles the group labels
(preserving group sizes), and recomputes the statistics from raw counts
using the same counting routine as the fit. The tested context set is
the union of the contexts the groups' trees actually represent. For two
groups the behavioral measure is JSD, which is **not** the symmetric-KL
distance used by
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md);
the `distance_matrix` component does use
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)
(symmetric KL) for consistency with the pairwise function.

## References

Lin, J. (1991). Divergence measures based on the Shannon entropy. *IEEE
Transactions on Information Theory*, 37(1), 145-151.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society: Series B*, 57(1), 289-300.

## See also

[`compare_trees`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for the pairwise permutation test.

## Examples

``` r
set.seed(1)
gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                           prob = c(.2,.6,.2)), simplify = FALSE)
gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                           prob = c(.2,.2,.6)), simplify = FALSE)
grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                    max_depth = 1L)
cmp <- compare_groups(grp, iter = 199L)
cmp
#> <transitiontrees_group_comparison>  2 groups, 199 permutations
#>   groups   : x, y
#>   contexts : 4 tested
#> 
#> omnibus (permutation):
#>        axis                 statistic  value p_value
#>  behavioral count-weighted JSD (bits) 184.09   0.005
#>       usage                   sum G^2 208.68   0.005
#> 
#> top pathways by behavioral divergence (JSD):
#>  pathway depth count_total modal_x modal_y flips jsd_bits jsd_padj usage_padj
#>        B     1         226       B       C  TRUE    0.171    0.005      0.007
#>  (start)     0         640       B       C  TRUE    0.162    0.005         NA
#>        A     1         117       B       C  TRUE    0.142    0.005      0.050
#>        C     1         217       B       C  TRUE    0.116    0.005      0.007
as.data.frame(cmp)
#>   pathway depth count_total count_x count_y modal_x modal_y flips  jsd_bits
#> 1       B     1         226     175      51       B       C  TRUE 0.1706917
#> 2 (start)     0         640     320     320       B       C  TRUE 0.1619964
#> 3       A     1         117      49      68       B       C  TRUE 0.1418159
#> 4       C     1         217      56     161       B       C  TRUE 0.1163298
#>   jsd_p jsd_padj   usage_g2 usage_p usage_padj
#> 1 0.005    0.005 119.098496   0.005     0.0075
#> 2 0.005    0.005         NA      NA         NA
#> 3 0.005    0.005   3.914325   0.050     0.0500
#> 4 0.005    0.005  85.666722   0.005     0.0075
```
