# Draw a forward trajectory tree

`plot_trajectories()` draws a forward trajectory tree of the sequences a
tree was fitted on: it starts from a common root and follows each
sequence state by state from its beginning, so every node is a prefix
(the first state, the first two states, and so on) and its position on
the x-axis is the move number. A context tree reads sequences backwards,
each node a history that ends at the moment of prediction; the
trajectory tree reads them forwards. It shows how sequences fan out from
their first state, which opening paths are common, and where they
diverge.

The tree is colored in one of two ways. With `measure = "frequency"` the
color of a node is the number of sequences that share the prefix. With
`measure = "predictability"` it is the probability the fitted model
gives to the node's last state after the preceding prefix,
\\P(\mathrm{state} \mid \mathrm{history})\\, with the history truncated
to the tree's maximum depth and read through
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md);
the root distribution is used for first states. Edge width is the number
of sequences in both cases, and darker color is a higher value. A pruned
tree shows predictability under the pruned model.

Only prefixes shared by at least `min_count` sequences are drawn. For
the backward context tree use
[`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md);
for the memory needed by one pathway use
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md).

## Usage

``` r
plot_trajectories(
  tree,
  measure = c("frequency", "predictability"),
  min_count = 4L
)
```

## Arguments

- tree:

  A `transitiontrees` object that carries its sequences.

- measure:

  `"frequency"` (default) or `"predictability"`.

- min_count:

  Smallest number of sequences sharing a prefix for it to be drawn.
  Default 4.

## Value

A ggplot object: the prefix tree, move number on the x-axis.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object, or carries no
  sequences.

- `transitiontrees_empty_result`:

  (error) no prefix is shared by `min_count` sequences.

## See also

[`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md)
for the backward context tree;
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md);
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
plot_trajectories(tree, min_count = 10)

plot_trajectories(prune_tree(tree), measure = "predictability",
                  min_count = 10)
```
