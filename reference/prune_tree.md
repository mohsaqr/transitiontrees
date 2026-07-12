# Prune a Context Tree

Removes nodes that do not earn their depth under the chosen criterion.
Pruning is applied bottom-up: a node is dropped when extending its
parent's prediction at this context produces less information /
likelihood than the depth penalty allows.

Renamed from `prune()` to avoid collision with other `prune()` generics.

## Usage

``` r
prune_tree(
  tree,
  criterion = c("G2", "KL", "AIC", "BIC"),
  alpha = 0.05,
  threshold = 0.005
)
```

## Arguments

- tree:

  A `transitiontrees`, or a `transitiontrees_group`, in which case each
  member is pruned and the group wrapper is preserved.

- criterion:

  One of `"G2"` (likelihood-ratio test against parent; default), `"KL"`
  (per-context Kullback-Leibler against parent), `"AIC"` (Akaike
  penalty), `"BIC"` (Bayesian penalty). Case-sensitive.

- alpha:

  Numeric in (0, 1). Significance level for `"G2"`; ignored otherwise.
  Default 0.05.

- threshold:

  Numeric. Minimum information gain in nats for `"KL"`; ignored
  otherwise. Default 0.005.

## Value

A pruned `transitiontrees` with `tree$pruned = TRUE` and `tree$pruning`
carrying the criterion + threshold settings.

## Details

For each leaf, compute the criterion against its parent. If the
criterion does not exceed its threshold, drop the leaf and revisit the
parent. Repeat until stable. The root is never dropped. Surviving nodes
keep their original smoothed `prob` vector (whatever smoothing scheme
was applied at fit time).

Note on units: the `"KL"` `threshold` is in **nats** (natural log),
whereas the `divergence` column reported by
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
/
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md)
is in **bits** (log base 2). Multiply a nats threshold by `1 / log(2)`
(~1.4427) to read it on the pathway-table scale.

## References

Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia. *Machine
Learning*, 25, 117-149.

## Examples

``` r
# \donttest{
seqs   <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                    simplify = FALSE)
tree   <- context_tree(seqs, max_depth = 4)
pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)
# }
```
