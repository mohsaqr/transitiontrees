# Student engagement state sequences (stslist)

A wide set of student engagement-state sequences as a `stslist` (the
state-sequence object `seqdef()` produces). Used to demonstrate loading
sequence objects directly into
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).
Bundled example dataset.

## Usage

``` r
engagement
```

## Format

A `stslist` with 1000 rows (learners) and 25 columns (time-steps).
States `"Active"`, `"Average"`, `"Disengaged"`; `"%"` marks missing/void
positions.

## Source

Bundled example dataset.

## Examples

``` r
data(engagement)
context_tree(engagement, max_depth = 2L)
#> <transitiontrees>  17 nodes, depth <= 2, 4 states  [unpruned]
#>   alphabet : %, Active, Average, Disengaged
#>   fit on   : 1000 sequences, 25000 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=25000  -> Active (0.50)
#> |-- %         n=145    -> % (1.00)
#> |   |-- Active    n=88     -> % (1.00)
#> |   |-- Average   n=15     -> % (1.00)
#> |   `-- Disengaged  n=42     -> % (1.00)
#> |-- Active    n=12063  -> Active (0.85)
#> |   |-- Active    n=9915   -> Active (0.85)
#> |   |-- Average   n=1490   -> Active (0.85)
#> |   `-- Disengaged  n=315    -> Active (0.84)
#> |-- Average   n=5118   -> Average (0.54)
#> |   |-- Active    n=1044   -> Average (0.53)
#> |   |-- Average   n=2699   -> Average (0.53)
#> |   `-- Disengaged  n=1056   -> Average (0.56)
#> `-- Disengaged  n=6674   -> Disengaged (0.78)
#>     |-- Active    n=586    -> Disengaged (0.81)
#>     |-- Average   n=723    -> Disengaged (0.73)
#>     `-- Disengaged  n=5027   -> Disengaged (0.78) 
```
