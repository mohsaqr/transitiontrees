# Compare Pruning Criteria on One Tree

Prunes a fitted tree under several criteria — holding `alpha` and
`threshold` fixed — and returns a tidy one-row-per-criterion summary of
how aggressively each trims the tree. A convenience wrapper over
repeated
[`prune_tree`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
calls that collapses the usual
[`vapply()`](https://rdrr.io/r/base/lapply.html) criterion loop into one
call.

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

  A `transitiontrees` (typically unpruned).

- criterion:

  Character vector of criteria to compare. Defaults to all four: `"G2"`,
  `"KL"`, `"AIC"`, `"BIC"`.

- alpha:

  Significance level for `"G2"` (and the AIC/BIC penalties' chi-square
  cutoff). Default 0.05.

- threshold:

  Minimum information gain in nats for `"KL"`. Default 0.005.

## Value

A `data.frame` with one row per criterion (in the order given by
`criterion`) and columns `criterion`, `n_nodes` (post-prune size) and
`reduction_pct` (percent of the original nodes removed).

## See also

[`prune_tree`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
to apply one criterion,
[`tune_tree`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
for cross-validated selection.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- replicate(80, sample(c("A", "B", "C"), 14, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 4L, min_count = 3L)
compare_pruning(tree)
#>   criterion n_nodes reduction_pct
#> 1        G2      10          91.7
#> 2        KL     118           2.5
#> 3       AIC      25          79.3
#> 4       BIC      24          80.2
compare_pruning(tree, criterion = c("G2", "BIC"), alpha = 0.01)
#>   criterion n_nodes reduction_pct
#> 1        G2       1          99.2
#> 2       BIC      24          80.2
# }
```
