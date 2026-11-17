# Plot pathways as a heatmap of next-state probabilities

`plot_pathways()` draws the pathway table of a context tree as a
heatmap: one row per pathway, one column per next state, and the fitted
probability \\P(\mathrm{next} \mid \mathrm{pathway})\\ as fill. The
pathway table lists, for every history, how often it occurs and what
tends to follow it. The most likely next state of each row is printed in
bold, other probabilities of at least 0.05 in plain text, and pathways
whose most likely next state differs from that of their shorter suffix
are marked with a leading `>`. One display then shows which pathways are
common, which lead almost always to one state, which leave the next
state open, and where a longer history changes the prediction that a
first-order model would make.

Rows are the `top` pathways by count, by divergence from the shorter
suffix, or by depth (the root is not shown). For the full distribution
of a few contexts as bar charts use
[`plot_distributions()`](https://pak.dynasite.org/transitiontrees/reference/plot_distributions.md);
to rank pathways by how much their history matters use
[`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md).

## Usage

``` r
plot_pathways(
  tree,
  top = 20L,
  sort_by = c("count", "divergence", "depth"),
  min_count = 5L,
  show_flips = TRUE,
  title = NULL,
  ...
)
```

## Arguments

- tree:

  A `transitiontrees` object.

- top:

  Maximum number of pathways shown. Default 20.

- sort_by:

  Order of the rows: `"count"` (default), `"divergence"` or `"depth"`
  (deepest first).

- min_count:

  Pathways with a lower count are dropped. Default 5.

- show_flips:

  Mark pathways whose most likely next state differs from their
  parent's. Default `TRUE`.

- title:

  Plot title. Default `NULL`: a title naming the ordering.

- ...:

  Ignored.

## Value

A ggplot object: a heatmap of pathways by next states.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object.

- `transitiontrees_empty_result`:

  (error) no pathway reaches `min_count`.

## See also

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
for the table;
[`plot_distributions()`](https://pak.dynasite.org/transitiontrees/reference/plot_distributions.md),
[`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
plot_pathways(tree, top = 12)

plot_pathways(tree, sort_by = "divergence", top = 12)
```
