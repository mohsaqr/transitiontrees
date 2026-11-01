# Summary of a Context Tree

Summary of a Context Tree

## Usage

``` r
# S3 method for class 'transitiontrees'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees`.

- ...:

  Ignored.

## Value

A `summary.transitiontrees` object. The `$table` slot is the canonical
pathway data.frame from
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
(columns `pathway`, `depth`, `count`, `likely_next`, `next_probability`,
`divergence`, `changes_prediction`), re-sorted by `(depth, -count)` so
the structural tree order is read top-to-bottom. Its other components
are the scalars `n_states`, `max_depth`, `n_nodes` and `pruned`.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2L)
summary(tree)
#> <transitiontrees summary>  13 nodes, depth <= 2, 3 states  [unpruned]
#> 
#>                   pathway depth count likely_next next_probability divergence
#>                   (start)     0  1870     Average        0.4347594         NA
#>                   Average     1   751     Average        0.6098535 0.11356246
#>                    Active     1   658      Active        0.6975684 0.34948716
#>                Disengaged     1   325  Disengaged        0.4830769 0.40306556
#>          Active -> Active     2   433      Active        0.7852194 0.02860157
#>        Average -> Average     2   419     Average        0.6778043 0.01466796
#>         Active -> Average     2   160     Average        0.5187500 0.12282588
#>         Average -> Active     2   144      Active        0.5000000 0.14727560
#>  Disengaged -> Disengaged     2   139  Disengaged        0.6762590 0.10977363
#>     Average -> Disengaged     2   134     Average        0.5000000 0.04611817
#>  changes_prediction
#>                  NA
#>               FALSE
#>                TRUE
#>                TRUE
#>               FALSE
#>               FALSE
#>               FALSE
#>               FALSE
#>               FALSE
#>                TRUE
#> # ... 3 more rows (use as.data.frame(tree) for the full table)
```
