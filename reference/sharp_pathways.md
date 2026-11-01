# Sharpest Pathways

Returns the top `n` pathways by predictive sharpness – the probability
mass on their modal next state. High values indicate strongly
deterministic continuations; low values indicate ambiguous next-state
distributions.

## Usage

``` r
sharp_pathways(tree, top = 10L, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees`.

- top:

  Integer. Number of pathways to return. Default 10.

- min_count:

  Integer. Drop pathways with fewer than this many occurrences. Default
  1.

## Value

A data.frame, same columns as
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
sorted by `next_probability` descending.

## Examples

``` r
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
sharp_pathways(tree, top = 5)
#>       pathway depth count likely_next next_probability divergence
#> 1 A -> B -> B     3    20           B        0.7000000 0.27863441
#> 2 A -> B -> A     3    22           A        0.6363636 0.07652960
#> 3 B -> C -> B     3     8           A        0.6241250 0.49292214
#> 4 A -> A -> B     3    20           B        0.6000000 0.15408490
#> 5 C -> A -> B     3    16           A        0.5625000 0.06185401
#>   changes_prediction
#> 1              FALSE
#> 2              FALSE
#> 3              FALSE
#> 4               TRUE
#> 5              FALSE
```
