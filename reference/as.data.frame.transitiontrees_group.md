# Convert a group of transition trees to one data frame of contexts

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
the context table of each tree in a `transitiontrees_group`, as returned
by
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
under a leading `group` column, so that the contexts of all groups can
be read, filtered or plotted together.

## Usage

``` r
# S3 method for class 'transitiontrees_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_group` object.

- row.names, optional:

  Not used.

- ...:

  Passed to
  [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md):
  `min_count`, `sort_by` and `decreasing`.

## Value

A `data.frame` with one row per group and context: a `group` column
followed by the columns of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md);
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
to test group differences context by context.

## Examples

``` r
groups <- context_tree(group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2)
head(as.data.frame(groups, min_count = 50))
#>   group           pathway depth count likely_next next_probability  divergence
#> 1  High           (start)     0 13721   consensus        0.2660885          NA
#> 2  High         consensus     1  3392        plan        0.3643868 0.341524645
#> 3  High              plan     1  2892        plan        0.3278008 0.196424995
#> 4  High           discuss     1  2003   consensus        0.4248627 0.459573925
#> 5  High           emotion     1  1550   consensus        0.3361290 0.509322703
#> 6  High consensus -> plan     2  1158        plan        0.3298791 0.001560704
#>   changes_prediction
#> 1                 NA
#> 2               TRUE
#> 3               TRUE
#> 4              FALSE
#> 5              FALSE
#> 6              FALSE
```
