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

A pruned `transitiontrees` (same structure as
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
returns, fewer nodes) whose `pruned` component is `TRUE` and whose
`pruning` component records `criterion`, `alpha` and `threshold`; the
print banner reports both. A `transitiontrees_group` input returns a
`transitiontrees_group` of pruned members.

## Details

For each leaf, compute the criterion against its parent. If the
criterion does not exceed its threshold, drop the leaf and revisit the
parent. Repeat until stable. The root is never dropped. `"G2"` keeps a
leaf when its likelihood-ratio statistic against the parent exceeds the
chi-square critical value at `alpha` with \\\|S\| - 1\\ degrees of
freedom (one test per leaf, no multiplicity correction). `"AIC"` and
`"BIC"` keep a leaf when its penalised log-likelihood beats the
parent's, with penalty \\2(\|S\| - 1)\\ and \\\log(n)(\|S\| - 1)\\
respectively (\\n\\ = the leaf's count). Surviving nodes keep their
original smoothed `prob` vector (whatever smoothing scheme was applied
at fit time).

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

Akaike, H. (1974). A new look at the statistical model identification.
*IEEE Transactions on Automatic Control*, 19(6), 716-723.

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics*, 6(2), 461-464.

## Examples

``` r
set.seed(1)
seqs   <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                    simplify = FALSE)
tree   <- context_tree(seqs, max_depth = 4)
pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)
pruned
#> <transitiontrees>  8 nodes, depth <= 4, 3 states  [pruned]
#>   alphabet : A, B, C
#>   fit on   : 50 sequences, 600 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#>   pruned by: G2   alpha = 0.05
#> (start)   n=600    -> A (0.37)
#> |-- A         n=200    -> A (0.38)
#> |   `-- A         n=66     -> A (0.38)
#> |       `-- B         n=31     -> C (0.39)
#> |           `-- B         n=7      -> A (0.71)
#> `-- B         n=181    -> A (0.40)
#>     `-- B         n=58     -> B (0.43)
#>         `-- A         n=20     -> B (0.70) 
prune_tree(tree, criterion = "BIC")
#> <transitiontrees>  26 nodes, depth <= 4, 3 states  [pruned]
#>   alphabet : A, B, C
#>   fit on   : 50 sequences, 600 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#>   pruned by: BIC   alpha = 0.05
#> (start)   n=600    -> A (0.37)
#> |-- A         n=200    -> A (0.38)
#> |   `-- A         n=66     -> A (0.38)
#> |       `-- B         n=31     -> C (0.39)
#> |           `-- B         n=7      -> A (0.71)
#> |-- B         n=181    -> A (0.40)
#> |   |-- A         n=56     -> A (0.43)
#> |   |   `-- B         n=15     -> A (0.33)
#> |   |       `-- C         n=5      -> A (0.60)
#> |   |-- B         n=58     -> B (0.43)
#> |   |   `-- A         n=20     -> B (0.70)
#> |   `-- C         n=51     -> A (0.43)
#> |       |-- A         n=18     -> A (0.39)
#> |       |   `-- C         n=5      -> A (0.80)
#> |       |-- B         n=8      -> A (0.62)
#> |       `-- C         n=18     -> A (0.39)
#> |           `-- C         n=5      -> B (0.60)
#> `-- C         n=169    -> C (0.35)
#>     |-- B         n=42     -> C (0.43)
#>     |   |-- B         n=14     -> A (0.43)
#>     |   |   `-- B         n=8      -> A (0.62)
#>     |   `-- C         n=14     -> C (0.50)
#>     |       `-- C         n=5      -> C (0.60)
#>     `-- C         n=51     -> B (0.37)
#>         `-- B         n=15     -> A (0.47)
#>             `-- C         n=5      -> A (0.80) 
```
