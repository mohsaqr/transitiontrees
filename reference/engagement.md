# Student engagement sequences (wide format)

`engagement` is a synthetic dataset of the engagement states of 1,000
students over a study program, generated to resemble the engagement
trajectories analyzed by Saqr, López-Pernas, Helske and Hrastinski
(2023). It has one row per student and one column per time point (23 to
25 observed states per student). Each state is `"Active"`, `"Average"`
or `"Disengaged"`; missing values occur only at the end of a row.

## Usage

``` r
engagement
```

## Format

A `data.frame` with 1,000 rows (students) and 25 character columns `T1`
to `T25`, holding `"Active"`, `"Average"`, `"Disengaged"` or `NA`.

## Source

The `engagement` sequence object of the tna package, converted to
character with its void and missing codes set to `NA` (identical
values).

## References

Saqr, M., López-Pernas, S., Helske, S., & Hrastinski, S. (2023). The
longitudinal association between engagement and achievement varies by
time, students' profiles, and achievement state: A full program study.
*Computers & Education*, 199, 104787.
[doi:10.1016/j.compedu.2023.104787](https://doi.org/10.1016/j.compedu.2023.104787)

## See also

[`trajectories`](https://pak.dynasite.org/transitiontrees/reference/trajectories.md)
for a smaller engagement dataset;
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Examples

``` r
tree <- context_tree(engagement, max_depth = 3, min_count = 20)
tree
#> <transitiontrees>  39 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 20
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |   |-- Active    n=7994   -> Active (0.86)
#> |   |   |-- Average   n=1215   -> Active (0.85)
#> |   |   `-- Disengaged  n=256    -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   |   |-- Active    n=296    -> Active (0.86)
#> |   |   |-- Average   n=781    -> Active (0.85)
#> |   |   `-- Disengaged  n=316    -> Active (0.87)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |       |-- Active    n=22     -> Active (0.73)
#> |       |-- Average   n=37     -> Active (0.84)
#> |       `-- Disengaged  n=238    -> Active (0.86)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |   |-- Active    n=859    -> Average (0.54)
#> |   |   |-- Average   n=118    -> Average (0.53)
#> |   |   `-- Disengaged  n=29     -> Average (0.59)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   |   |-- Active    n=520    -> Average (0.54)
#> |   |   |-- Average   n=1394   -> Average (0.53)
#> |   |   `-- Disengaged  n=567    -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> |       |-- Active    n=79     -> Average (0.66)
#> |       |-- Average   n=140    -> Average (0.59) 
#> ... 13 more nodes (use as.data.frame(x) or summary(x))
divergent_pathways(tree, top = 5)
#>                           pathway depth count likely_next next_probability
#> 1                      Disengaged     1  6598  Disengaged        0.7859958
#> 2                          Active     1 11894      Active        0.8590886
#> 3                         Average     1  5063     Average        0.5478965
#> 4  Active -> Disengaged -> Active     3    22      Active        0.7272727
#> 5 Disengaged -> Active -> Average     3    29     Average        0.5862069
#>   divergence changes_prediction
#> 1 0.94531802               TRUE
#> 2 0.41831831              FALSE
#> 3 0.38363675               TRUE
#> 4 0.25459907              FALSE
#> 5 0.05596631              FALSE
```
