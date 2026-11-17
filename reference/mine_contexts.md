# Find the contexts that make a given next state likely or unlikely

`mine_contexts()` finds the contexts that make a given next state likely
or unlikely. It scans every context of a fitted tree for the probability
it gives one chosen next state and returns the contexts whose
probability lies in a requested range, from the highest probability
down. This serves questions about one outcome state, such as after which
histories disengagement is most likely, or after which histories a
learner almost never returns to planning.

Each row reports the context, its depth and count, the probability of
the chosen state after it, and whether that state is the context's most
likely next state. Probabilities in rare contexts are imprecise; set
`min_count` to keep only contexts with enough occurrences. To compare
the state's probability across groups, fit a grouped tree and use
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

## Usage

``` r
mine_contexts(tree, state, min_prob = NULL, max_prob = NULL, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees` object.

- state:

  The next state of interest, one of the tree's states.

- min_prob, max_prob:

  Keep contexts whose probability of `state` is at least `min_prob` and
  at most `max_prob`. `NULL` (default) leaves that bound open.

- min_count:

  Keep contexts that occurred at least this many times. Default 1.

## Value

A `data.frame` with one row per context and columns `pathway`, `depth`,
`count`, `state`, `prob` (probability of `state` after the context) and
`is_modal` (whether `state` is the most likely next state), sorted by
`prob` and then `count`, both decreasing. No matching context gives a
zero-row table with these columns.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or `state`
is not a single one of its states.

## See also

[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
[`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md),
[`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
mine_contexts(tree, state = "Disengaged", min_prob = 0.4, min_count = 10)
#>                                  pathway depth count      state      prob
#> 1 Disengaged -> Disengaged -> Disengaged     3    81 Disengaged 0.7037037
#> 2               Disengaged -> Disengaged     2   139 Disengaged 0.6762590
#> 3    Average -> Disengaged -> Disengaged     3    45 Disengaged 0.6000000
#> 4       Average -> Disengaged -> Average     3    64 Disengaged 0.4843750
#> 5                             Disengaged     1   325 Disengaged 0.4830769
#> 6       Average -> Average -> Disengaged     3    62 Disengaged 0.4032258
#> 7        Disengaged -> Active -> Average     3    10 Disengaged 0.4000000
#>   is_modal
#> 1     TRUE
#> 2     TRUE
#> 3     TRUE
#> 4     TRUE
#> 5     TRUE
#> 6    FALSE
#> 7    FALSE
mine_contexts(tree, state = "Active", max_prob = 0.05)
#>                               pathway depth count  state       prob is_modal
#> 1 Average -> Disengaged -> Disengaged     3    45 Active 0.02222222    FALSE
#> 2 Disengaged -> Disengaged -> Average     3    31 Active 0.00100000    FALSE
#> 3     Active -> Disengaged -> Average     3     7 Active 0.00100000    FALSE
```
