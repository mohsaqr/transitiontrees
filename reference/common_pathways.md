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
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
common_pathways(tree, top = 8)
#>   pathway depth count likely_next next_probability  divergence
#> 1 (start)     0   600           A        0.3666667          NA
#> 2       A     1   200           A        0.3800000 0.003825681
#> 3       B     1   181           A        0.3977901 0.013246404
#> 4       C     1   169           C        0.3491124 0.008133498
#> 5  B -> A     2    67           A        0.4776119 0.030761672
#> 6  A -> A     2    66           A        0.3787879 0.003552395
#> 7  A -> C     2    60           B        0.3500000 0.001283401
#> 8  B -> B     2    58           B        0.4310345 0.021097272
#>   changes_prediction
#> 1                 NA
#> 2              FALSE
#> 3              FALSE
#> 4               TRUE
#> 5              FALSE
#> 6              FALSE
#> 7               TRUE
#> 8               TRUE
common_pathways(tree, top = 8, depth = 3L)   # restrict to depth-3
#>        pathway depth count likely_next next_probability divergence
#> 14 B -> A -> A     3    31           C        0.3870968 0.02360426
#> 15 A -> B -> A     3    22           A        0.6363636 0.07652960
#> 16 B -> B -> B     3    22           A        0.4090909 0.13625381
#> 17 C -> B -> A     3    21           A        0.3809524 0.03090333
#> 18 A -> A -> B     3    20           B        0.6000000 0.15408490
#> 19 A -> B -> B     3    20           B        0.7000000 0.27863441
#> 20 B -> A -> C     3    19           C        0.4736842 0.06159616
#> 21 C -> C -> B     3    18           A        0.3888889 0.01229457
#>    changes_prediction
#> 14               TRUE
#> 15              FALSE
#> 16               TRUE
#> 17              FALSE
#> 18               TRUE
#> 19              FALSE
#> 20               TRUE
#> 21              FALSE
```
