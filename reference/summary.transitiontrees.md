# Summarize a fitted transition tree

[`summary()`](https://rdrr.io/r/base/summary.html) summarizes a fitted
tree by its size (numbers of states and contexts, depth, pruned or not)
and its table of contexts ordered by depth and, within depth, by
decreasing count, so that the table reads from the root downward.

## Usage

``` r
# S3 method for class 'transitiontrees'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees` object.

- ...:

  Not used.

## Value

An object of class `summary.transitiontrees`, which prints a header and
the first rows of the context table (columns as in
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)).
For the full table use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
tree; for fit statistics use
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md).

## See also

[`print.summary.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/print.summary.transitiontrees.md),
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
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
