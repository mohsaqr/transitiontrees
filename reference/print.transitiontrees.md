# Print a fitted transition tree

[`print()`](https://rdrr.io/r/base/print.html) prints a header (numbers
of contexts, depth and states, the data the tree was fitted on, the
smoothing scheme and `min_count`, and the pruning criterion when the
tree is pruned), followed by the tree drawn as text. Each line of the
drawing is one context: the state it adds to its parent's history (the
older state), its count `n`, and the most likely next state with its
probability. Long trees are cut after `max_lines` lines;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives
every context.

## Usage

``` r
# S3 method for class 'transitiontrees'
print(x, max_lines = 25L, digits = 2L, ...)
```

## Arguments

- x:

  A `transitiontrees` object.

- max_lines:

  Maximum number of tree lines to print. Default 25.

- digits:

  Decimal places for probabilities. Default 2.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees.md),
[`as.data.frame.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees.md),
[`plot.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
print(tree, max_lines = 10)
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
