# Count how often each state occurs

`state_frequencies()` tabulates the states of the sequences a tree was
fitted on. The counts are those of the root context, which conditions on
no history and counts every observed position once; they form the
marginal distribution of states, against which every longer context is
compared. For a `transitiontrees_group`, the table is computed for each
group.

## Usage

``` r
state_frequencies(tree)
```

## Arguments

- tree:

  A `transitiontrees` or a `transitiontrees_group` object.

## Value

A `data.frame` with one row per state (per group and state, with a
leading `group` column, for a `transitiontrees_group`) and columns
`state`, `count` (observed positions, weighted for a weighted fit) and
`proportion`, sorted by decreasing count, with ties in the order of the
tree's states.

## Conditions

`transitiontrees_bad_input` when `tree` is neither a fitted tree nor a
group of them.

## See also

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
state_frequencies(tree)
#>        state count proportion
#> 1    Average   813  0.4347594
#> 2     Active   703  0.3759358
#> 3 Disengaged   354  0.1893048
```
