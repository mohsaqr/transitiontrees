# Mine Contexts by Next-State Probability

Scan every context in a fitted tree for a chosen next state and return
those whose predicted probability for that state falls in a requested
range. A tidy context-mining table: "in which histories is the next move
`state` unusually likely or unlikely?"

## Usage

``` r
mine_contexts(tree, state, min_prob = NULL, max_prob = NULL, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees`.

- state:

  Character. The next state to score, one of the tree's alphabet.

- min_prob, max_prob:

  Numeric in \\\[0, 1\]\\ or `NULL`. Keep contexts whose
  `P(state | context)` is at least `min_prob` and/or at most `max_prob`.
  `NULL` (default) leaves that side unbounded.

- min_count:

  Integer. Drop contexts with fewer than this many occurrences. Default
  1.

## Value

A data.frame with columns `pathway`, `depth`, `count`, `state`, `prob`
(`P(state | context)`), and `is_modal` (whether `state` is the context's
most likely next state; ties broken by alphabet order), sorted by `prob`
descending. The empty case returns a 0-row data.frame with the same
schema.

## Examples

``` r
set.seed(1)
seqs <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 2L)
mine_contexts(tree, state = "A", min_prob = 0.4)
#>   pathway depth count state      prob is_modal
#> 1  A -> B     2    53     A 0.4528302     TRUE
#> 2  B -> A     2    61     A 0.4426230     TRUE
#> 3  C -> B     2    52     A 0.4230769     TRUE
```
