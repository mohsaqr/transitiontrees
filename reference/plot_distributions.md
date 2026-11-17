# Plot the next-state distribution of each context

`plot_distributions()` draws the full next-state distribution
\\P(\mathrm{next} \mid \mathrm{context})\\ of selected contexts as small
multiples, one bar chart per context with the most likely next state in
a darker shade. The display shows, for each history, whether the next
state is close to determined (one tall bar) or open (bars of similar
height), and which alternatives remain. By default the most frequent
contexts are shown.

The same probabilities are shown for many contexts at once as a heatmap
by
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
and along one history's chain of shorter suffixes by
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md).

## Usage

``` r
plot_distributions(tree, contexts = NULL, top = 12L, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees` object.

- contexts:

  Contexts to show, in arrow notation. Default `NULL`: the `top` most
  frequent contexts.

- top:

  Number of contexts when `contexts` is `NULL`. Default 12.

- min_count:

  Contexts with a lower count are dropped. Default 1.

## Value

A ggplot object with one bar-chart panel per context.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object.

- `transitiontrees_empty_result`:

  (error) no context is left to plot.

## See also

[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md);
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
for the table of contexts.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
plot_distributions(tree, top = 6)

plot_distributions(tree, contexts = c("Active", "Active -> Disengaged"))
```
