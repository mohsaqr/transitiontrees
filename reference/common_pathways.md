# List the most frequent pathways of a tree

`common_pathways()` returns the pathways that occurred most often,
optionally at one depth. Frequent pathways are the typical routes
through the process: the histories most learners or sessions actually
pass through, and the state each of them most often leads to. Frequency
says nothing about whether a pathway's prediction differs from that of
its shorter suffix; for that, see
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md).

## Usage

``` r
common_pathways(tree, top = 10L, depth = NULL, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees` object.

- top:

  Number of pathways to return. Default 10.

- depth:

  Return only pathways of this depth. `NULL` (default) returns all
  depths.

- min_count:

  Keep pathways that occurred at least this many times. Default 1.

## Value

A `data.frame` with the columns of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
one row per pathway, sorted by `count` in decreasing order, at most
`top` rows.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md),
[`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md),
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
common_pathways(tree, top = 5)
#>              pathway depth count likely_next next_probability divergence
#> 1            (start)     0  1870     Average        0.4347594         NA
#> 2            Average     1   751     Average        0.6098535 0.11356246
#> 3             Active     1   658      Active        0.6975684 0.34948716
#> 4   Active -> Active     2   433      Active        0.7852194 0.02860157
#> 5 Average -> Average     2   419     Average        0.6778043 0.01466796
#>   changes_prediction
#> 1                 NA
#> 2              FALSE
#> 3               TRUE
#> 4              FALSE
#> 5              FALSE
common_pathways(tree, top = 5, depth = 3)
#>                                  pathway depth count likely_next
#> 1             Active -> Active -> Active     3   316      Active
#> 2          Average -> Average -> Average     3   248     Average
#> 3 Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#> 4           Active -> Average -> Average     3    80     Average
#> 5            Active -> Active -> Average     3    74     Average
#>   next_probability  divergence changes_prediction
#> 1        0.8354430 0.011491875              FALSE
#> 2        0.7177419 0.006242261              FALSE
#> 3        0.7037037 0.008106021              FALSE
#> 4        0.5750000 0.101876320              FALSE
#> 5        0.5540541 0.024836699              FALSE
```
