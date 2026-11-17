# List the pathways with the most certain next state

`sharp_pathways()` returns the pathways whose most likely next state has
the highest probability. A pathway with a next-state probability close
to 1 is followed almost always by the same state, a near-deterministic
step in the process; a probability close to \\1/k\\ (\\k\\ states) means
the next state is hard to predict from that history. Sharp pathways with
small counts are often artefacts of few observations; set `min_count`
accordingly.

## Usage

``` r
sharp_pathways(tree, top = 10L, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees` object.

- top:

  Number of pathways to return. Default 10.

- min_count:

  Keep pathways that occurred at least this many times. Default 1.

## Value

A `data.frame` with the columns of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
one row per pathway, sorted by `next_probability` in decreasing order,
at most `top` rows.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md)
for the entropy of each prediction.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
sharp_pathways(tree, top = 5, min_count = 20)
#>                                  pathway depth count likely_next
#> 1             Active -> Active -> Active     3   316      Active
#> 2                       Active -> Active     2   433      Active
#> 3          Average -> Average -> Average     3   248     Average
#> 4    Disengaged -> Disengaged -> Average     3    31     Average
#> 5 Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#>   next_probability  divergence changes_prediction
#> 1        0.8354430 0.011491875              FALSE
#> 2        0.7852194 0.028601567              FALSE
#> 3        0.7177419 0.006242261              FALSE
#> 4        0.7085484 0.198404690              FALSE
#> 5        0.7037037 0.008106021              FALSE
```
