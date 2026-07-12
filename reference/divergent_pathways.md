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
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
divergent_pathways(tree, top = 6)
#>       pathway depth count likely_next next_probability divergence
#> 1 A -> A -> B     3    12           A        0.5825833  0.5737596
#> 2 C -> A -> B     3    19           C        0.5263158  0.2240645
#> 3 A -> B -> A     3    14           B        0.5714286  0.2173227
#> 4 B -> C -> A     3    17           A        0.5882353  0.1732741
#> 5 B -> B -> A     3    17           A        0.4705882  0.1606417
#> 6 C -> C -> B     3    12           B        0.5000000  0.1107270
#>   changes_prediction
#> 1               TRUE
#> 2               TRUE
#> 3               TRUE
#> 4              FALSE
#> 5              FALSE
#> 6              FALSE
divergent_pathways(tree, flips_only = TRUE)
#>        pathway depth count likely_next next_probability divergence
#> 1  A -> A -> B     3    12           A        0.5825833 0.57375959
#> 2  C -> A -> B     3    19           C        0.5263158 0.22406455
#> 3  A -> B -> A     3    14           B        0.5714286 0.21732269
#> 7  C -> A -> A     3    21           A        0.4761905 0.10217770
#> 8  A -> C -> A     3    21           B        0.4285714 0.09413230
#> 11 B -> A -> C     3    13           B        0.5384615 0.06615899
#> 12 C -> B -> A     3    13           C        0.4615385 0.06590266
#> 15 C -> B -> C     3    18           A        0.4444444 0.05257821
#> 17 A -> A -> A     3    20           B        0.4000000 0.04762184
#> 18 B -> C -> B     3    23           C        0.4347826 0.04672169
#>    changes_prediction
#> 1                TRUE
#> 2                TRUE
#> 3                TRUE
#> 7                TRUE
#> 8                TRUE
#> 11               TRUE
#> 12               TRUE
#> 15               TRUE
#> 17               TRUE
#> 18               TRUE
# }
```
