# Most Common Pathways in a Fitted Tree

Returns the top `n` pathways by occurrence count – the trajectories the
data actually contains many copies of.

## Usage

``` r
common_pathways(tree, top = 10L, depth = NULL, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees`.

- top:

  Integer. Number of pathways to return. Default 10.

- depth:

  Integer or NULL. Restrict to pathways of this exact depth. `NULL`
  (default) keeps all depths.

- min_count:

  Integer. Minimum count cut-off. Default 1.

## Value

A data.frame, same columns as
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
sorted by count descending.

## Examples

``` r
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
common_pathways(tree, top = 8)
#>   pathway depth count likely_next next_probability   divergence
#> 1 (start)     0   600           B        0.3500000           NA
#> 2       A     1   193           A        0.3523316 0.0002749574
#> 3       B     1   190           A        0.3578947 0.0007928251
#> 4       C     1   167           B        0.4071856 0.0101555958
#> 5  A -> A     2    65           C        0.3538462 0.0094907810
#> 6  B -> A     2    64           A        0.3593750 0.0012072641
#> 7  A -> B     2    62           A        0.4193548 0.0146403185
#> 8  B -> B     2    61           B        0.3442623 0.0044815150
#>   changes_prediction
#> 1                 NA
#> 2               TRUE
#> 3               TRUE
#> 4              FALSE
#> 5               TRUE
#> 6              FALSE
#> 7              FALSE
#> 8               TRUE
common_pathways(tree, top = 8, depth = 3L)   # restrict to depth-3
#>        pathway depth count likely_next next_probability  divergence
#> 14 A -> B -> A     3    25           A        0.4000000 0.005868202
#> 15 B -> A -> A     3    23           C        0.4782609 0.089534432
#> 16 A -> C -> A     3    20           A        0.4000000 0.026560781
#> 17 B -> A -> B     3    20           A        0.4000000 0.013307423
#> 18 C -> B -> A     3    20           B        0.3500000 0.011654583
#> 19 B -> B -> B     3    20           C        0.5000000 0.091211911
#> 20 A -> A -> A     3    20           C        0.4000000 0.018680728
#> 21 A -> A -> C     3    19           A        0.4210526 0.007066973
#>    changes_prediction
#> 14              FALSE
#> 15              FALSE
#> 16              FALSE
#> 17              FALSE
#> 18               TRUE
#> 19               TRUE
#> 20              FALSE
#> 21              FALSE
# }
```
