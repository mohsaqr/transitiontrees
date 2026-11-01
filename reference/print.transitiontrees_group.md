# Print a Group of Context Trees

Print a Group of Context Trees

## Usage

``` r
# S3 method for class 'transitiontrees_group'
print(x, ...)
```

## Arguments

- x:

  A `transitiontrees_group` (named list of `transitiontrees`s).

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
grp <- context_tree(trajectories, group = rep(c("early", "late"), each = 69),
                    max_depth = 2L)
print(grp)
#> <transitiontrees_group>  2 groups
#>   early           13 nodes, depth <= 2, 69 seq, 969 obs
#>   late            13 nodes, depth <= 2, 67 seq, 901 obs
```
