# Measure how much each step of history changes the prediction

`tree_dependence()` measures how much each step of history changes the
prediction of a fitted tree, that is, where the process is path
dependent. A process is path dependent where what happens next depends
on more than the current state. In a transition tree this shows as a
context whose next-state distribution differs from that of its parent,
the same context without its oldest state. The difference has two sides:
how far the distribution moves (its divergence) and whether it becomes
more or less certain (its change in entropy).

It returns one row for every context other than the root, with the
Kullback-Leibler divergence of its next-state distribution from its
parent's, the Shannon entropy of both distributions and their
difference, and the most likely next state under each (Cover & Thomas,
2006). A positive `entropy_drop` means the additional state of history
makes the next state more predictable; a negative one means it makes it
less predictable. The divergence is the quantity that
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
compares with its threshold under `criterion = "KL"` (there in nats).

For the pathway table with counts and ranked views, use
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
and
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md);
to see the full next-state distributions along one pathway, use
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md).

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

  A `transitiontrees` object.

- base:

  Base of the logarithm: 2 (default) for bits, `exp(1)` for nats.

- sort_by:

  Column to sort by, in decreasing order: `"divergence"` (default),
  `"entropy_drop"`, `"entropy"`, `"count"` or `"depth"`.

- top:

  Number of rows to keep after sorting. `NULL` (default) keeps all.

## Value

A `data.frame` with one row per non-root context whose parent is in the
tree, and columns `pathway`, `depth`, `count`, `divergence` (from the
parent's distribution), `entropy` (of this context's distribution),
`entropy_before` (of the parent's), `entropy_drop` (`entropy_before`
minus `entropy`), `likely_next`, `likely_before` (the parent's most
likely next state) and `changes_prediction`. A tree with only a root
returns a zero-row table with these columns.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## References

Cover, T. M., & Thomas, J. A. (2006). *Elements of information theory*
(2nd ed.). Wiley.

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md),
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
[`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
tree_dependence(tree, top = 5)
#>                           pathway depth count divergence   entropy
#> 1  Active -> Disengaged -> Active     3     6  0.6773610 0.6623258
#> 2 Disengaged -> Active -> Average     3    10  0.5477701 1.3609640
#> 3                      Disengaged     1   325  0.4030656 1.4032668
#> 4 Disengaged -> Average -> Active     3    12  0.3933342 0.6623258
#> 5                          Active     1   658  0.3494872 1.0404534
#>   entropy_before entropy_drop likely_next likely_before changes_prediction
#> 1       1.361016  0.698689876      Active        Active              FALSE
#> 2       1.365499  0.004535442     Average       Average              FALSE
#> 3       1.507629  0.104361741  Disengaged       Average               TRUE
#> 4       1.124941  0.462615346     Average        Active               TRUE
#> 5       1.507629  0.467175189      Active       Average               TRUE
tree_dependence(prune_tree(tree), sort_by = "entropy_drop")
#>                                pathway depth count divergence   entropy
#> 1  Disengaged -> Disengaged -> Average     3    31 0.19840469 0.8802163
#> 2                               Active     1   658 0.34948716 1.0404534
#> 3      Disengaged -> Average -> Active     3    12 0.39333418 0.6623258
#> 4             Disengaged -> Disengaged     2   139 0.10977363 1.1558077
#> 5                     Active -> Active     2   433 0.02860157 0.8490387
#> 6                              Average     1   751 0.11356246 1.3544863
#> 7                   Average -> Average     2   419 0.01466796 1.2289362
#> 8                           Disengaged     1   325 0.40306556 1.4032668
#> 9      Disengaged -> Active -> Average     3    10 0.54777014 1.3609640
#> 10               Disengaged -> Average     2   122 0.17618640 1.3516073
#> 11                   Active -> Average     2   160 0.12282588 1.3654995
#> 12               Average -> Disengaged     2   134 0.04611817 1.4301429
#> 13                   Average -> Active     2   144 0.14727560 1.1249411
#> 14        Active -> Average -> Average     3    80 0.10187632 1.3380583
#> 15     Active -> Average -> Disengaged     3    17 0.31329586 1.5485652
#> 16                Active -> Disengaged     2    23 0.34780951 1.5743318
#> 17         Average -> Active -> Active     3    70 0.23300616 1.2048347
#>    entropy_before entropy_drop likely_next likely_before changes_prediction
#> 1       1.3516073  0.471390982     Average       Average              FALSE
#> 2       1.5076285  0.467175189      Active       Average               TRUE
#> 3       1.1249411  0.462615346     Average        Active               TRUE
#> 4       1.4032668  0.247459099  Disengaged    Disengaged              FALSE
#> 5       1.0404534  0.191414680      Active        Active              FALSE
#> 6       1.5076285  0.153142262     Average       Average              FALSE
#> 7       1.3544863  0.125550085     Average       Average              FALSE
#> 8       1.5076285  0.104361741  Disengaged       Average               TRUE
#> 9       1.3654995  0.004535442     Average       Average              FALSE
#> 10      1.3544863  0.002878951     Average       Average              FALSE
#> 11      1.3544863 -0.011013211     Average       Average              FALSE
#> 12      1.4032668 -0.026876133     Average    Disengaged               TRUE
#> 13      1.0404534 -0.084487794      Active        Active              FALSE
#> 14      1.2289362 -0.109122155     Average       Average              FALSE
#> 15      1.4301429 -0.118422293      Active       Average               TRUE
#> 16      1.4032668 -0.171064986      Active    Disengaged               TRUE
#> 17      0.8490387 -0.355796060      Active        Active              FALSE
```
