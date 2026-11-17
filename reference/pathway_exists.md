# Test whether a pathway is a context of the tree

`pathway_exists()` reports whether a history is itself a context (node)
of the fitted tree, that is, whether it occurred at least `min_count`
times and, for a pruned tree, survived pruning. A history that is not a
context is still predicted, from its longest suffix in the tree; see
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md).

## Usage

``` r
pathway_exists(tree, pathway)
```

## Arguments

- tree:

  A `transitiontrees` object.

- pathway:

  The history, oldest state first: a string in arrow notation or a
  character vector.

## Value

A single logical value.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## See also

[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
[`subtree()`](https://pak.dynasite.org/transitiontrees/reference/subtree.md),
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
pathway_exists(tree, "Active -> Disengaged -> Active")
#> [1] TRUE
pathway_exists(prune_tree(tree), "Active -> Disengaged -> Active")
#> [1] FALSE
```
