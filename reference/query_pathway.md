# Look up the next-state distribution after a pathway

`query_pathway()` returns the probability distribution that the fitted
tree assigns to the next state after a given history, or the probability
of one particular next state. For example, it gives the probability of
Disengaged after the history Active, Average, read directly from the
fitted model.

By default the history need not be a context of the tree: as in
prediction, the longest suffix of the history that is in the tree is
used, falling back to the root when no suffix matches. With
`exact = TRUE` only a context of the tree is accepted, and other
histories return `NA`. To find which context matched, use
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md)
with `tidy = TRUE`; to test whether a pathway is a context, use
[`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md).

## Usage

``` r
query_pathway(tree, pathway, next_state = NULL, exact = FALSE)
```

## Arguments

- tree:

  A `transitiontrees` object.

- pathway:

  The history, oldest state first: a single string in arrow notation
  (`"A -> B"`) or a character vector (`c("A", "B")`). `"(start)"` or an
  empty vector denotes the root.

- next_state:

  A state whose probability is returned. `NULL` (default) returns the
  whole distribution.

- exact:

  Require `pathway` to be a context of the tree. Default `FALSE` (use
  the longest matching suffix).

## Value

With `next_state`, a single probability. Otherwise a named numeric
vector with one probability per state, in the order of the tree's
states. `NA` for a pathway that is not a context when `exact = TRUE`.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`next_state` is not one of its states.

## See also

[`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md),
[`subtree()`](https://pak.dynasite.org/transitiontrees/reference/subtree.md),
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md),
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
query_pathway(tree, "Active -> Average")
#>     Active    Average Disengaged 
#>    0.37500    0.51875    0.10625 
query_pathway(tree, c("Active", "Average"), next_state = "Disengaged")
#> [1] 0.10625
```
