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
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
sharp_pathways(tree, top = 5)
#>       pathway depth count likely_next next_probability divergence
#> 1 C -> A -> C     3    11           B        0.6363636 0.13121800
#> 2 C -> B -> B     3    15           A        0.6000000 0.13746158
#> 3 C -> A -> A     3    17           A        0.5882353 0.12401427
#> 4 B -> C -> A     3    19           B        0.5789474 0.04083662
#> 5 A -> C -> B     3    21           C        0.5714286 0.15093086
#>   changes_prediction
#> 1              FALSE
#> 2              FALSE
#> 3              FALSE
#> 4              FALSE
#> 5              FALSE
# }
```
