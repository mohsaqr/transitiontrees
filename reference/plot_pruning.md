# Plot how much memory one pathway needs

`plot_pruning()` plots how much memory one pathway needs: how many past
states are needed to predict what follows it. A context tree keeps a
longer history only when it predicts the next state differently from its
shorter suffix, and the plot makes that decision visible for one
pathway. It follows the pathway's suffix chain (the full history, then
the history without its earliest state, and so on down to the root) and
draws the next-state distribution at every link as its own panel,
deepest history on the left and the root on the right.

Each context is classified by the \\G^2\\ rule of
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md):
the likelihood-ratio statistic of its next-state counts against its
one-shorter parent is compared with the chi-square critical value at
`alpha`. A context whose own \\G^2\\ exceeds it *adds memory* and is
drawn solid; a context below it that is kept only because a deeper
context exceeds it is drawn at medium opacity; a context below it with
no deeper context exceeding it is *pruned* and drawn faded. Each panel
title gives the context, the decision and \\G^2\\. The distributions are
the node values reported by
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
and
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md).

Pass an unpruned tree so that the whole chain is visible. To prune the
whole tree use
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md);
to compare pruning rules use
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md).

## Usage

``` r
plot_pruning(tree, pathway, alpha = 0.05)
```

## Arguments

- tree:

  A `transitiontrees` object, usually unpruned.

- pathway:

  One fitted context in arrow notation (`"A -> B -> C"`, earliest state
  on the left).

- alpha:

  Significance level of the \\G^2\\ test. Default 0.05.

## Value

A ggplot object with one panel per context of the suffix chain.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object; `alpha` is not in
  (0, 1); `pathway` is not a context of the tree.

## See also

[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md),
[`plot_distributions()`](https://pak.dynasite.org/transitiontrees/reference/plot_distributions.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
plot_pruning(tree, "Average -> Active -> Active")
```
