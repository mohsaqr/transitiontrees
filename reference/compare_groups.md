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

  A `transitiontrees_group`.

- iter:

  Integer. Number of label permutations. Default 999.

- min_count:

  Integer. Drop contexts whose total count across all groups is below
  this. Default 1.

- seed:

  Integer or `NULL`. RNG seed. Default 1.

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

  Two-row data.frame: the behavioral and usage global statistics with
  permutation p-values.

- distance_matrix:

  K x K symmetric-KL distance matrix between the groups (from
  [`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)).

- groups, iter, seed, n_contexts:

  Configuration.

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

## See also

[`compare_trees`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for the pairwise permutation test.

## Examples

``` r
# \donttest{
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
#>        axis                 statistic   value p_value
#>  behavioral count-weighted JSD (bits) 206.318   0.005
#>       usage                   sum G^2 235.718   0.005
#> 
#> top pathways by behavioral divergence (JSD):
#>  pathway depth count_total modal_x modal_y flips jsd_bits jsd_padj usage_padj
#>        C     1         213       B       C  TRUE    0.217    0.005      0.007
#>  (start)     0         640       B       C  TRUE    0.189    0.005         NA
#>        B     1         224       B       C  TRUE    0.129    0.005      0.007
#>        A     1         123       B       C  TRUE    0.085    0.005      1.000
# }
```
