# Coerce a Group of Trees to One Tidy Data Frame

Row-binds each group's
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
table, tagged with a leading `group` column, so the whole batch is one
tidy frame ready to filter, facet, or join.

## Usage

``` r
# S3 method for class 'transitiontrees_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_group`.

- row.names, optional:

  Ignored.

- ...:

  Forwarded to
  [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Value

A data.frame with one row per (group, pathway): a leading `group` column
identifying the source tree, then the canonical pathway columns of
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Examples

``` r
grp <- context_tree(trajectories, group = rep(c("early", "late"), each = 69),
                    max_depth = 2L)
as.data.frame(grp)
#>    group                  pathway depth count likely_next next_probability
#> 1  early                  (start)     0   969     Average        0.4592363
#> 2  early                  Average     1   412     Average        0.6383495
#> 3  early                   Active     1   290      Active        0.6793103
#> 4  early       Average -> Average     2   240     Average        0.7000000
#> 5  early               Disengaged     1   198  Disengaged        0.5050505
#> 6  early         Active -> Active     2   183      Active        0.7814208
#> 7  early Disengaged -> Disengaged     2    91  Disengaged        0.6703297
#> 8  early    Average -> Disengaged     2    80     Average        0.5000000
#> 9  early    Disengaged -> Average     2    71     Average        0.5492958
#> 10 early        Active -> Average     2    70     Average        0.5285714
#> 11 early        Average -> Active     2    61     Average        0.5081967
#> 12 early     Disengaged -> Active     2    21      Active        0.5238095
#> 13 early     Active -> Disengaged     2    15  Disengaged        0.4000000
#> 14  late                  (start)     0   901      Active        0.4328524
#> 15  late                   Active     1   368      Active        0.7119565
#> 16  late                  Average     1   339     Average        0.5752212
#> 17  late         Active -> Active     2   250      Active        0.7880000
#> 18  late       Average -> Average     2   179     Average        0.6480447
#> 19  late               Disengaged     1   127  Disengaged        0.4488189
#> 20  late        Active -> Average     2    90     Average        0.5111111
#> 21  late        Average -> Active     2    83      Active        0.5301205
#> 22  late    Average -> Disengaged     2    54     Average        0.5000000
#> 23  late    Disengaged -> Average     2    51     Average        0.4705882
#> 24  late Disengaged -> Disengaged     2    48  Disengaged        0.6875000
#> 25  late     Active -> Disengaged     2     8      Active        0.5000000
#> 26  late     Disengaged -> Active     2     8      Active        0.4995000
#>    divergence changes_prediction
#> 1          NA                 NA
#> 2  0.11860117              FALSE
#> 3  0.41368908               TRUE
#> 4  0.01235085              FALSE
#> 5  0.32508480               TRUE
#> 6  0.03863114              FALSE
#> 7  0.08068567              FALSE
#> 8  0.06324216               TRUE
#> 9  0.17816331              FALSE
#> 10 0.23939796              FALSE
#> 11 0.18535737               TRUE
#> 12 0.11611845              FALSE
#> 13 0.20764597              FALSE
#> 14         NA                 NA
#> 15 0.28445668              FALSE
#> 16 0.10462308               TRUE
#> 17 0.02161718              FALSE
#> 18 0.01625662              FALSE
#> 19 0.50429558               TRUE
#> 20 0.04536041              FALSE
#> 21 0.12209663              FALSE
#> 22 0.03028335               TRUE
#> 23 0.17987959              FALSE
#> 24 0.17403527              FALSE
#> 25 0.77249882               TRUE
#> 26 0.19341184              FALSE
```
