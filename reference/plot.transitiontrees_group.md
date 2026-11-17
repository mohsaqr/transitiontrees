# Plot the context tree of each group

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws every
tree of a `transitiontrees_group` with
[`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md),
one figure per group with the group name as caption, so that the memory
structure of the groups can be compared side by side. Each figure is
printed as it is made. To test or map the differences between groups use
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
and
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md).

## Usage

``` r
# S3 method for class 'transitiontrees_group'
plot(x, ...)
```

## Arguments

- x:

  A `transitiontrees_group` object.

- ...:

  Passed to
  [`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md),
  such as `style`.

## Value

Invisibly, a named list with one plot per group, in group order: ggplot
objects, or htmlwidgets for `style = "interactive"`.

## See also

[`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md);
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
and
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for the differences between groups.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
plot(grp, show_prediction = FALSE)

```
