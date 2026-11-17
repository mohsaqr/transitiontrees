# Extract the branch of a tree below one context

`subtree()` returns the part of a fitted tree that lies below a given
context: the context itself and every longer context that ends with it,
that is, every history that extends it further into the past. The branch
shows how the prediction after a given recent history is refined by what
happened before it, and can be summarized, converted to a data frame and
plotted like a full tree.

Contexts keep their full names, and the states, smoothing and other
settings are those of the original tree. The printed header names the
context at which the branch was cut.

## Usage

``` r
subtree(tree, pathway)
```

## Arguments

- tree:

  A `transitiontrees` object.

- pathway:

  The context at which to cut, oldest state first: a string in arrow
  notation or a character vector. It must be a context of the tree.

## Value

A `transitiontrees` object holding `pathway` and its descendants.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`pathway` is not one of its contexts.

## See also

[`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md),
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
[`plot.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
branch <- subtree(tree, "Disengaged")
summary(branch)
#> <transitiontrees summary>  11 nodes, depth <= 3, 3 states  [unpruned]
#> 
#>                                 pathway depth count likely_next
#>                              Disengaged     1   325  Disengaged
#>                Disengaged -> Disengaged     2   139  Disengaged
#>                   Average -> Disengaged     2   134     Average
#>                    Active -> Disengaged     2    23      Active
#>  Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#>        Average -> Average -> Disengaged     3    62     Average
#>     Disengaged -> Average -> Disengaged     3    46     Average
#>     Average -> Disengaged -> Disengaged     3    45  Disengaged
#>         Active -> Average -> Disengaged     3    17      Active
#>          Active -> Active -> Disengaged     3     9      Active
#>  next_probability  divergence changes_prediction
#>         0.4830769          NA                 NA
#>         0.6762590 0.109773627              FALSE
#>         0.5000000 0.046118168               TRUE
#>         0.3913043 0.347809514               TRUE
#>         0.7037037 0.008106021              FALSE
#>         0.5161290 0.026852413              FALSE
#>         0.6086957 0.048844466              FALSE
#>         0.6000000 0.079870032              FALSE
#>         0.4117647 0.313295862               TRUE
#>         0.4444444 0.024571720              FALSE
#> # ... 1 more rows (use as.data.frame(tree) for the full table)
as.data.frame(branch, sort_by = "depth")
#>                                   pathway depth count likely_next
#> 1          Active -> Active -> Disengaged     3     9      Active
#> 2      Active -> Disengaged -> Disengaged     3     7  Disengaged
#> 3  Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#> 4     Average -> Disengaged -> Disengaged     3    45  Disengaged
#> 5        Average -> Average -> Disengaged     3    62     Average
#> 6     Disengaged -> Average -> Disengaged     3    46     Average
#> 7         Active -> Average -> Disengaged     3    17      Active
#> 8                    Active -> Disengaged     2    23      Active
#> 9                Disengaged -> Disengaged     2   139  Disengaged
#> 10                  Average -> Disengaged     2   134     Average
#> 11                             Disengaged     1   325  Disengaged
#>    next_probability  divergence changes_prediction
#> 1         0.4444444 0.024571720              FALSE
#> 2         0.5714286 0.312801206              FALSE
#> 3         0.7037037 0.008106021              FALSE
#> 4         0.6000000 0.079870032              FALSE
#> 5         0.5161290 0.026852413              FALSE
#> 6         0.6086957 0.048844466              FALSE
#> 7         0.4117647 0.313295862               TRUE
#> 8         0.3913043 0.347809514               TRUE
#> 9         0.6762590 0.109773627              FALSE
#> 10        0.5000000 0.046118168               TRUE
#> 11        0.4830769          NA                 NA
```
