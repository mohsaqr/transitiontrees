# Count the contexts of a transition tree

`n_nodes()` returns the number of contexts (nodes) in a fitted tree, the
root included: the number shown in the printed header and a direct
measure of the tree's complexity, for example before and after pruning.

## Usage

``` r
n_nodes(tree)
```

## Arguments

- tree:

  A `transitiontrees` or `transitiontrees_group` object.

## Value

A single integer; for a `transitiontrees_group`, a named integer vector
with one count per group.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or group.

## See also

[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md),
[`summary.transitiontrees_group()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
n_nodes(tree)
#> [1] 38
n_nodes(prune_tree(tree))
#> [1] 18
```
