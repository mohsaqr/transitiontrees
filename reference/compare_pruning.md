# Compare pruning criteria on one tree

`compare_pruning()` prunes one fitted tree under each of several
criteria and reports how many contexts each keeps. The table shows how
strongly the criteria differ on the data at hand before one of them is
applied with
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md);
the criteria themselves are described there. A smaller tree is not by
itself a better one: to compare pruned trees by how well they predict
new sequences, use
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

## Usage

``` r
compare_pruning(
  tree,
  criterion = c("G2", "KL", "AIC", "BIC"),
  alpha = 0.05,
  threshold = 0.005
)
```

## Arguments

- tree:

  A `transitiontrees` object, usually unpruned.

- criterion:

  Criteria to compare. Default all four: `"G2"`, `"KL"`, `"AIC"`,
  `"BIC"`.

- alpha:

  Significance level for `"G2"`. Default 0.05.

- threshold:

  Minimum divergence in nats for `"KL"`. Default 0.005.

## Value

A `data.frame` with one row per criterion, in the order of `criterion`,
and columns `criterion`, `n_nodes` (contexts kept, root included) and
`reduction_pct` (percentage of the original contexts removed).

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree,
`criterion` is not a non-empty character vector, or `alpha` or
`threshold` is invalid.

## See also

[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md),
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`compare_smoothing()`](https://pak.dynasite.org/transitiontrees/reference/compare_smoothing.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 4)
compare_pruning(tree)
#>   criterion n_nodes reduction_pct
#> 1        G2      25          69.5
#> 2        KL      78           4.9
#> 3       AIC      42          48.8
#> 4       BIC      36          56.1
compare_pruning(tree, criterion = c("G2", "BIC"), alpha = 0.01)
#>   criterion n_nodes reduction_pct
#> 1        G2      21          74.4
#> 2       BIC      36          56.1
```
