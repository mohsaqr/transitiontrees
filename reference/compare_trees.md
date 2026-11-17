# Test whether two context trees differ, with a permutation test

`compare_trees()` tests whether two fitted context trees, one per group,
predict differently overall, with a permutation test. The test statistic
is the count-weighted symmetric Kullback-Leibler divergence between the
trees' next-state distributions, the quantity returned by
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md);
two groups of sequences can pass through the same states and still
differ in what tends to follow a given history.

The divergence has no known sampling distribution and is positive even
when both groups come from one process, so the reference distribution is
built by permutation. The sequences of the two trees are pooled,
reassigned at random to two groups of the original sizes, and two trees
are refitted with each original tree's depth, minimum count, smoothing
and pruning; this is repeated `iter` times. The p-value is the share of
permuted distances at least as large as the observed one, counting the
observed value once (\\(1 + b) / (iter + 1)\\). A per-context table
shows which histories contribute most to the distance. The result is a
`transitiontrees_comparison` object, read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) and drawn with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

The comparison is pairwise and treats sequences as exchangeable. For
more than two groups, for separate tests of what follows a context and
how often it is reached, for per-context p-values with false discovery
rate control, or for a permutation stratified by person when one person
contributes several sequences, use
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

## Usage

``` r
compare_trees(tree_a, tree_b = NULL, iter = 200L, seed = 1L, symmetric = TRUE)
```

## Arguments

- tree_a, tree_b:

  Unweighted `transitiontrees` objects fitted to the two groups.
  Alternatively a `transitiontrees_group` with two groups as `tree_a`
  and `tree_b = NULL` (default); its two trees are compared in group
  order.

- iter:

  Number of permutations. Default 200.

- seed:

  Seed for the permutations. Default 1. The call sets the session's
  random-number state.

- symmetric:

  `TRUE` (default) for the symmetric divergence; `FALSE` for \\D\_{KL}(A
  \\ B)\\ alone.

## Value

An object of class `transitiontrees_comparison`. Its tables are returned
by `as.data.frame(x, what = "pathways")` (one row per context in either
tree: counts in each tree and the two directed and the symmetric
divergence) and `as.data.frame(x, what = "null")` (one row per
permutation); [`summary()`](https://rdrr.io/r/base/summary.html) gives
one row with the observed distance, the null mean, standard deviation
and 95th percentile, and the p-value;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the null
distribution.

## Conditions

- `transitiontrees_bad_input`:

  (error) an argument is not a `transitiontrees` object; `tree_b`
  missing; a `transitiontrees_group` without exactly two groups, or
  given together with `tree_b`; `iter` not a positive number.

- `transitiontrees_weighted_unsupported`:

  (error) a tree was fitted with weights; the permutation refits from
  the raw sequences and would drop them.

- `transitiontrees_alphabet_mismatch`:

  (error) the trees have different alphabets.

## See also

[`summary.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_comparison.md),
[`as.data.frame.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_comparison.md)
and
[`plot.transitiontrees_comparison`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_comparison.md)
for the result;
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)
for the distance alone;
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for several groups and per-context tests.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
cmp <- compare_trees(grp, iter = 20)
cmp
#> <transitiontrees_comparison>  iter = 20
#>   observed distance : 0.0471
#>   null mean         : 0.00294
#>   p-value           : 0.0476
#> 
#> top divergent pathways:
#>     pathway count_a count_b divergence_ab divergence_ba divergence_sym
#>   synthesis     278     374         0.300         0.311          0.305
#>     discuss    2003    1948         0.178         0.204          0.191
#>    cohesion     938     757         0.120         0.115          0.117
#>       adapt     141     368         0.044         0.094          0.069
#>   consensus    3392    2937         0.050         0.046          0.048
#>  coregulate     894    1076         0.046         0.047          0.046
summary(cmp)
#> <transitiontrees_comparison summary>
#>  distance null_mean  null_sd null_q95 p_value iter n_pathways
#>    0.0471   0.00294 0.000556  0.00378  0.0476   20         10
```
