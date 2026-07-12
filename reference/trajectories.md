# Student engagement trajectories

An example set of categorical learning-engagement trajectories used in
the transitiontrees examples and the “trajectories” vignette. Each row
is one learner, each column a time-step, and each cell the engagement
state at that step. Trailing `NA`s mark the end of a trajectory. This
wide character matrix is exactly the shape
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
consumes.

## Usage

``` r
trajectories
```

## Format

A character matrix with 138 rows (learners) and 15 columns (time-steps).
Three states: `"Active"`, `"Average"`, `"Disengaged"`.

## Source

Bundled example dataset.

## Examples

``` r
data(trajectories)
dim(trajectories)
#> [1] 138  15
tree <- context_tree(trajectories)
tree
#> <transitiontrees>  145 nodes, depth <= 5, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |   |-- Active    n=240    -> Active (0.87)
#> |   |   |   |   |-- Active    n=192    -> Active (0.87)
#> |   |   |   |   |-- Average   n=17     -> Active (0.76)
#> |   |   |   |   `-- Disengaged  n=7      -> Active (0.71)
#> |   |   |   |-- Average   n=37     -> Active (0.59)
#> |   |   |   |   |-- Active    n=12     -> Active (0.83)
#> |   |   |   |   `-- Average   n=18     -> Active (0.50)
#> |   |   |   `-- Disengaged  n=10     -> Active (0.90)
#> |   |   |       `-- Active    n=5      -> Active (0.80)
#> |   |   |-- Average   n=70     -> Active (0.53)
#> |   |   |   |-- Active    n=22     -> Active (0.55)
#> |   |   |   |   |-- Active    n=11     -> Active (0.64)
#> |   |   |   |   `-- Average   n=10     -> Average (0.50)
#> |   |   |   `-- Average   n=37     -> Active (0.49)
#> |   |   |       |-- Active    n=15     -> Active (0.53)
#> |   |   |       |-- Average   n=15     -> Active (0.47)
#> |   |   |       `-- Disengaged  n=7      -> Average (0.57)
#> |   |   `-- Disengaged  n=15     -> Active (0.67)
#> |   |       |-- Active    n=5      -> Active (1.00)
#> |   |       `-- Average   n=7      -> Average (0.57)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   |   |-- Active    n=53     -> Average (0.55) 
#> ... 119 more nodes (use as.data.frame(x) or summary(x))
```
