# Student engagement state sequences (wide format)

A wide set of student engagement-state sequences, one learner per row
and one time-step per column. Used to demonstrate fitting from wide
sequence data and from the network objects (tna / Nestimate) built on
top of it. Bundled example dataset.

## Usage

``` r
engagement
```

## Format

A `data.frame` with 1000 rows (learners) and 25 columns (time-steps).
States `"Active"`, `"Average"`, `"Disengaged"`; `NA` marks missing /
dropped-out positions.

## Source

Bundled example dataset.

## Examples

``` r
data(engagement)
context_tree(engagement, max_depth = 2L)
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> `-- Disengaged  n=6598   -> Disengaged (0.79)
#>     |-- Active    n=581    -> Disengaged (0.82)
#>     |-- Average   n=719    -> Disengaged (0.73)
#>     `-- Disengaged  n=4960   -> Disengaged (0.79) 
```
