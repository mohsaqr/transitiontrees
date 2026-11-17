# Plot the context-level differences between groups

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
result in one of two ways. The default lollipop chart ranks the contexts
by the Jensen-Shannon divergence between the groups' next-state
distributions, with point size for the context's count and color
separating contexts whose adjusted p-value is below `alpha`; it shows
where, given the same history, the groups continue differently.
`style = "matrix"` draws the pairwise distances between groups as a
labeled heatmap.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
plot(x, style = c("divergence", "matrix"), top = 15L, alpha = 0.05, ...)
```

## Arguments

- x:

  A `transitiontrees_group_comparison` object.

- style:

  `"divergence"` (default) for the ranked contexts, or `"matrix"` for
  the between-group distance heatmap.

- top:

  Number of contexts shown by `"divergence"`. Default 15.

- alpha:

  False discovery rate cutoff for highlighting. Default 0.05.

- ...:

  Ignored.

## Value

A ggplot object: for `"divergence"`, one row per context with its
divergence in bits; for `"matrix"`, one tile per pair of groups with its
symmetric Kullback-Leibler distance printed.

## Conditions

- `transitiontrees_empty_result`:

  (error) for `"divergence"`, no context has a positive divergence.

## See also

[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md);
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
for the next-state differences within each context.

## Examples

``` r
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 1)
cmp <- compare_groups(grp, iter = 20)
plot(cmp, alpha = 0.1)

plot(cmp, style = "matrix")
```
