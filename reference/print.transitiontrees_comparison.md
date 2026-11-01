# Print a context tree Comparison

Print a context tree Comparison

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
print(x, digits = 3L, n = 6L, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison` object.

- digits:

  Integer. Numeric digits for the printed table. Default 3.

- n:

  Integer. Number of top divergent pathways to print. Default 6.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
grp <- context_tree(trajectories, group = rep(c("early", "late"), each = 69),
                    max_depth = 2L)
print(compare_trees(grp, iter = 50L), n = 3L)
#> <transitiontrees_comparison>  iter = 50
#>   observed distance : 0.0246
#>   null mean         : 0.0173
#>   p-value           : 0.157
#> 
#> top divergent pathways:
#>                pathway count_a count_b divergence_ab divergence_ba
#>   Disengaged -> Active      21       8         0.599         0.173
#>   Active -> Disengaged      15       8         0.239         0.185
#>  Disengaged -> Average      71      51         0.050         0.067
#>  divergence_sym
#>           0.386
#>           0.212
#>           0.058
```
