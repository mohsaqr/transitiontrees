# Print a Context Tree

Print a Context Tree

## Usage

``` r
# S3 method for class 'transitiontrees'
print(x, max_lines = 25L, digits = 2L, ...)
```

## Arguments

- x:

  A `transitiontrees`.

- max_lines:

  Integer. Maximum tree-rendering lines. Default 25.

- digits:

  Integer. Probability digits. Default 2.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2L)
print(tree, max_lines = 10L)
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   `-- Disengaged  n=29     -> Active (0.52)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> `-- Disengaged  n=325    -> Disengaged (0.48)
#>     |-- Active    n=23     -> Active (0.39) 
#> ... 2 more nodes (use as.data.frame(x) or summary(x))
```
