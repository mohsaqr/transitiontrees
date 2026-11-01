# Coerce a context tree to a Tidy Data Frame

Returns the canonical tidy node table, identical to
`tree_pathways(tree)`: one row per context (pathway).

## Usage

``` r
# S3 method for class 'transitiontrees'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees`.

- row.names, optional:

  Ignored.

- ...:

  Forwarded to
  [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Value

A data.frame with columns `pathway`, `depth`, `count`, `likely_next`,
`next_probability`, `divergence`, `changes_prediction`; one row per
pathway. See
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
for the column meanings.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2L)
as.data.frame(tree)
#>                     pathway depth count likely_next next_probability divergence
#> 1                   (start)     0  1870     Average        0.4347594         NA
#> 2                   Average     1   751     Average        0.6098535 0.11356246
#> 3                    Active     1   658      Active        0.6975684 0.34948716
#> 4          Active -> Active     2   433      Active        0.7852194 0.02860157
#> 5        Average -> Average     2   419     Average        0.6778043 0.01466796
#> 6                Disengaged     1   325  Disengaged        0.4830769 0.40306556
#> 7         Active -> Average     2   160     Average        0.5187500 0.12282588
#> 8         Average -> Active     2   144      Active        0.5000000 0.14727560
#> 9  Disengaged -> Disengaged     2   139  Disengaged        0.6762590 0.10977363
#> 10    Average -> Disengaged     2   134     Average        0.5000000 0.04611817
#> 11    Disengaged -> Average     2   122     Average        0.5163934 0.17618640
#> 12     Disengaged -> Active     2    29      Active        0.5172414 0.12990434
#> 13     Active -> Disengaged     2    23      Active        0.3913043 0.34780951
#>    changes_prediction
#> 1                  NA
#> 2               FALSE
#> 3                TRUE
#> 4               FALSE
#> 5               FALSE
#> 6                TRUE
#> 7               FALSE
#> 8               FALSE
#> 9               FALSE
#> 10               TRUE
#> 11              FALSE
#> 12              FALSE
#> 13               TRUE
as.data.frame(tree, sort_by = "divergence")
#>                     pathway depth count likely_next next_probability divergence
#> 1                Disengaged     1   325  Disengaged        0.4830769 0.40306556
#> 2                    Active     1   658      Active        0.6975684 0.34948716
#> 3      Active -> Disengaged     2    23      Active        0.3913043 0.34780951
#> 4     Disengaged -> Average     2   122     Average        0.5163934 0.17618640
#> 5         Average -> Active     2   144      Active        0.5000000 0.14727560
#> 6      Disengaged -> Active     2    29      Active        0.5172414 0.12990434
#> 7         Active -> Average     2   160     Average        0.5187500 0.12282588
#> 8                   Average     1   751     Average        0.6098535 0.11356246
#> 9  Disengaged -> Disengaged     2   139  Disengaged        0.6762590 0.10977363
#> 10    Average -> Disengaged     2   134     Average        0.5000000 0.04611817
#> 11         Active -> Active     2   433      Active        0.7852194 0.02860157
#> 12       Average -> Average     2   419     Average        0.6778043 0.01466796
#> 13                  (start)     0  1870     Average        0.4347594         NA
#>    changes_prediction
#> 1                TRUE
#> 2                TRUE
#> 3                TRUE
#> 4               FALSE
#> 5               FALSE
#> 6               FALSE
#> 7               FALSE
#> 8               FALSE
#> 9               FALSE
#> 10               TRUE
#> 11              FALSE
#> 12              FALSE
#> 13                 NA
```
