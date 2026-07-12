# Per-Context Path Dependence of a Context Tree

For each non-root node in the tree, reports the Kullback-Leibler
divergence of its conditional next-state distribution against its
parent's. Large values flag contexts where extending memory by one more
step changes the prediction; the `changes_prediction` column flags
contexts where the most likely next state changes between the node and
its parent.

Renamed from `path_dependence()` to avoid a naming collision with a
sibling package.

## Usage

``` r
tree_dependence(
  tree,
  base = 2,
  sort_by = c("divergence", "entropy_drop", "entropy", "count", "depth"),
  top = NULL
)
```

## Arguments

- tree:

  A `transitiontrees`.

- base:

  Numeric. Logarithm base for the KL divergence. Default 2 (bits). Use
  `exp(1)` for nats or `10` for hartleys.

- sort_by:

  Character. Column to sort by, descending. One of `"divergence"`
  (default), `"entropy_drop"`, `"entropy"`, `"count"`, `"depth"`.

- top:

  Integer or `NULL`. If given, keep only the top `top` rows after
  sorting. Default `NULL` (all rows).

## Value

A data.frame with one row per non-root pathway, sorted by `divergence`
descending. Columns: `pathway`, `depth`, `count`, `divergence`
(Kullback-Leibler divergence from the parent's prediction), `entropy`
(Shannon entropy of this pathway's next-state distribution),
`entropy_before` (entropy of the parent's distribution), `entropy_drop`
(`entropy_before - entropy`, the uncertainty this step of history
removes), `likely_next` (this node's most likely next state),
`likely_before` (the parent context's most likely next state),
`changes_prediction` (`likely_next != likely_before`). The empty case
returns a 0-row data.frame with the same schema.

## Details

This is the *diagnostic* that the tree's pruning rule (under
`criterion = "KL"`) is comparing against its threshold. It answers the
substantive question: for which contexts does this tree disagree with a
memoryless / shorter-memory model, and where does that disagreement
actually flip the prediction?

The mean of `n * KL` across rows recovers, up to constants, the
chain-level mutual-information gain from the variable-depth model over
the order-1 model.

## References

Cover, T.M. & Thomas, J.A. (2006). *Elements of Information Theory*, 2nd
ed. Wiley.

## Examples

``` r
# \donttest{
seqs   <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                    simplify = FALSE)
tree   <- context_tree(seqs, max_depth = 3)
pruned <- prune_tree(tree, criterion = "G2")
tree_dependence(pruned)
#>       pathway depth count  divergence  entropy entropy_before entropy_drop
#> 1 C -> B -> C     3    11 0.415780328 1.004468       1.560781  0.556313425
#> 2      B -> C     2    44 0.019366728 1.560781       1.584554  0.023772726
#> 3           C     1   157 0.005768993 1.584554       1.577429 -0.007125206
#>   likely_next likely_before changes_prediction
#> 1           C             A               TRUE
#> 2           A             A              FALSE
#> 3           A             A              FALSE
# }
```
