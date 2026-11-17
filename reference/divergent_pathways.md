# List the pathways whose history most changes the prediction

`divergent_pathways()` returns the pathways whose next-state
distribution departs most from that of their parent, the same pathway
without its oldest state, measured by the Kullback-Leibler divergence in
bits. These are the places where remembering one more state of history
changes what is expected next, and therefore where the process depends
on more than its most recent state. With `flips_only = TRUE`, only
pathways whose most likely next state differs from the parent's are
returned.

Divergence is computed from smoothed probabilities and can be large for
rare pathways; set `min_count` to exclude them, and use
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
to check whether a divergence is stable under resampling. The root,
which has no parent, is never returned.

## Usage

``` r
divergent_pathways(tree, top = 10L, min_count = 1L, flips_only = FALSE)
```

## Arguments

- tree:

  A `transitiontrees` object.

- top:

  Number of pathways to return. Default 10.

- min_count:

  Keep pathways that occurred at least this many times. Default 1.

- flips_only:

  Return only pathways whose most likely next state differs from the
  parent's. Default `FALSE`.

## Value

A `data.frame` with the columns of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
one row per pathway, sorted by `divergence` in decreasing order, at most
`top` rows.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md),
[`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md),
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
divergent_pathways(tree, top = 5, min_count = 20)
#>                               pathway depth count likely_next next_probability
#> 1                          Disengaged     1   325  Disengaged        0.4830769
#> 2                              Active     1   658      Active        0.6975684
#> 3                Active -> Disengaged     2    23      Active        0.3913043
#> 4         Average -> Active -> Active     3    70      Active        0.5285714
#> 5 Disengaged -> Disengaged -> Average     3    31     Average        0.7085484
#>   divergence changes_prediction
#> 1  0.4030656               TRUE
#> 2  0.3494872               TRUE
#> 3  0.3478095               TRUE
#> 4  0.2330062              FALSE
#> 5  0.1984047              FALSE
divergent_pathways(tree, flips_only = TRUE)
#>                               pathway depth count likely_next next_probability
#> 3                          Disengaged     1   325  Disengaged        0.4830769
#> 4     Disengaged -> Average -> Active     3    12     Average        0.8318333
#> 5                              Active     1   658      Active        0.6975684
#> 6                Active -> Disengaged     2    23      Active        0.3913043
#> 7     Active -> Average -> Disengaged     3    17      Active        0.4117647
#> 10    Active -> Disengaged -> Average     3     7  Disengaged        0.5707143
#> 18 Disengaged -> Disengaged -> Active     3     9     Average        0.5555556
#> 23              Average -> Disengaged     2   134     Average        0.5000000
#> 24   Average -> Disengaged -> Average     3    64  Disengaged        0.4843750
#> 30       Average -> Active -> Average     3    64      Active        0.4531250
#>    divergence changes_prediction
#> 3  0.40306556               TRUE
#> 4  0.39333418               TRUE
#> 5  0.34948716               TRUE
#> 6  0.34780951               TRUE
#> 7  0.31329586               TRUE
#> 10 0.20137567               TRUE
#> 18 0.10602662               TRUE
#> 23 0.04611817               TRUE
#> 24 0.02940545               TRUE
#> 30 0.02076716               TRUE
```
