# Print a Context-Tree Summary

Print method for the object returned by
[`summary.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees.md):
a one-line banner followed by the leading rows of the canonical pathway
table.

## Usage

``` r
# S3 method for class 'summary.transitiontrees'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `summary.transitiontrees` object.

- n:

  Integer. Number of pathway rows to print. Default 10.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2L)
print(summary(tree), n = 5L)
#> <transitiontrees summary>  13 nodes, depth <= 2, 3 states  [unpruned]
#> 
#>           pathway depth count likely_next next_probability divergence
#>           (start)     0  1870     Average        0.4347594         NA
#>           Average     1   751     Average        0.6098535 0.11356246
#>            Active     1   658      Active        0.6975684 0.34948716
#>        Disengaged     1   325  Disengaged        0.4830769 0.40306556
#>  Active -> Active     2   433      Active        0.7852194 0.02860157
#>  changes_prediction
#>                  NA
#>               FALSE
#>                TRUE
#>                TRUE
#>               FALSE
#> # ... 8 more rows (use as.data.frame(tree) for the full table)
```
