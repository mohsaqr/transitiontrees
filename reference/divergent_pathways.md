# Most Predictively Divergent Pathways

Returns the top `n` pathways by Kullback-Leibler divergence from their
(k-1)-suffix. These are the pathways whose extended history adds the
most predictive information over the shorter one. Pathways whose most
likely next state actually flips between orders are marked in the
`changes_prediction` column.

## Usage

``` r
divergent_pathways(tree, top = 10L, min_count = 1L, flips_only = FALSE)
```

## Arguments

- tree:

  A `transitiontrees`.

- top:

  Integer. Number of pathways to return. Default 10.

- min_count:

  Integer. Drop pathways with fewer than this many occurrences. Default
  1.

- flips_only:

  Logical. If `TRUE`, return only pathways that flip the most likely
  next state between orders. Default `FALSE`.

## Value

A data.frame, same columns as
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
sorted by `divergence` descending.

## Examples

``` r
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
divergent_pathways(tree, top = 6)
#>       pathway depth count likely_next next_probability divergence
#> 1 B -> C -> B     3     8           A        0.6241250 0.49292214
#> 2 A -> B -> B     3    20           B        0.7000000 0.27863441
#> 3 A -> A -> B     3    20           B        0.6000000 0.15408490
#> 4 B -> B -> B     3    22           A        0.4090909 0.13625381
#> 5 C -> A -> C     3    17           A        0.4705882 0.11128385
#> 6 A -> B -> C     3     9           C        0.4444444 0.08368479
#>   changes_prediction
#> 1              FALSE
#> 2              FALSE
#> 3               TRUE
#> 4               TRUE
#> 5               TRUE
#> 6              FALSE
divergent_pathways(tree, flips_only = TRUE)
#>        pathway depth count likely_next next_probability divergence
#> 3  A -> A -> B     3    20           B        0.6000000 0.15408490
#> 4  B -> B -> B     3    22           A        0.4090909 0.13625381
#> 5  C -> A -> C     3    17           A        0.4705882 0.11128385
#> 9  B -> B -> C     3    14           A        0.4285714 0.06417903
#> 11 B -> A -> C     3    19           C        0.4736842 0.06159616
#> 12      C -> A     2    49           C        0.4081633 0.05968289
#> 17 B -> C -> C     3    15           A        0.4666667 0.04319042
#> 18 C -> C -> A     3    17           B        0.4117647 0.04042774
#> 26 B -> A -> A     3    31           C        0.3870968 0.02360426
#> 28      B -> B     2    58           B        0.4310345 0.02109727
#>    changes_prediction
#> 3                TRUE
#> 4                TRUE
#> 5                TRUE
#> 9                TRUE
#> 11               TRUE
#> 12               TRUE
#> 17               TRUE
#> 18               TRUE
#> 26               TRUE
#> 28               TRUE
```
