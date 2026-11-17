# Student engagement trajectories

`trajectories` is a coded dataset of the engagement states of 138
students at 15 successive time points, one row per student and one
column per time point. Each state is `"Active"`, `"Average"` or
`"Disengaged"`. Missing values occur only as a run at the end of a row
and mark that the student's trajectory ended (drop-out); two rows are
missing throughout and are ignored when a tree is fitted. The matrix has
the wide layout that
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
fits directly, and is small enough for fast examples.

## Usage

``` r
trajectories
```

## Format

A character matrix with 138 rows (students) and 15 columns (time points,
named `"1"` to `"15"`), holding `"Active"`, `"Average"`, `"Disengaged"`
or `NA`. 118 rows are complete; 1870 states are observed in all.

## Source

Bundled example dataset. The same matrix is distributed with the
Nestimate package as `trajectories`.

## See also

[`engagement`](https://pak.dynasite.org/transitiontrees/reference/engagement.md)
for a larger engagement dataset;
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
tree
#> <transitiontrees>  38 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |-- Average   n=70     -> Active (0.53)
#> |   |   `-- Disengaged  n=15     -> Active (0.67)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   |   |-- Active    n=53     -> Average (0.55)
#> |   |   |-- Average   n=66     -> Active (0.58)
#> |   |   `-- Disengaged  n=12     -> Average (0.83)
#> |   `-- Disengaged  n=29     -> Active (0.52)
#> |       |-- Active    n=6      -> Active (0.83)
#> |       |-- Average   n=13     -> Active (0.54)
#> |       `-- Disengaged  n=9      -> Average (0.56)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   |-- Active    n=74     -> Average (0.55)
#> |   |   |-- Average   n=64     -> Active (0.45)
#> |   |   `-- Disengaged  n=10     -> Average (0.50)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.57)
#> |   |   |-- Average   n=248    -> Average (0.72)
#> |   |   `-- Disengaged  n=63     -> Average (0.65)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> |       |-- Active    n=7      -> Disengaged (0.57)
#> |       |-- Average   n=64     -> Disengaged (0.48) 
#> ... 12 more nodes (use as.data.frame(x) or summary(x))
state_frequencies(tree)
#>        state count proportion
#> 1    Average   813  0.4347594
#> 2     Active   703  0.3759358
#> 3 Disengaged   354  0.1893048
```
