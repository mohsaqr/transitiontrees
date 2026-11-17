# Find the least and most typical sequences

`mine_sequences()` ranks sequences by how well a fitted tree predicts
them and returns the top `n`. With `which = "surprising"`, the sequences
with the highest perplexity come first: those whose transitions the tree
finds least expected, which are candidates for atypical cases or for
patterns that the fitting data did not contain. With
`which = "expected"`, the sequences with the lowest perplexity come
first: the most typical sequences under the tree.

The ranking is by per-sequence perplexity from
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
which adjusts for sequence length. When `newdata` are the fitting
sequences themselves, the ranking describes typicality within the data;
with other sequences it describes how well the tree generalizes to them.
To find where within a sequence the unexpected moves occur, use
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md).

## Usage

``` r
mine_sequences(tree, newdata, n = 10L, which = c("surprising", "expected"))
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Sequences to rank, in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- n:

  Number of sequences to return. Default 10.

- which:

  `"surprising"` (default; highest perplexity first) or `"expected"`
  (lowest perplexity first).

## Value

A `data.frame` with the columns of
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md)
(`sequence_id`, `n_scored`, `log_lik`, `perplexity`), one row per
sequence, at most `n` rows.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`newdata` has no usable sequence.

## See also

[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md),
[`mine_contexts()`](https://pak.dynasite.org/transitiontrees/reference/mine_contexts.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
mine_sequences(tree, engagement, n = 5)
#>   sequence_id n_scored   log_lik perplexity
#> 1         230       23 -30.83375   3.821328
#> 2         126       24 -31.93443   3.783317
#> 3         552       25 -32.46789   3.664587
#> 4         775       25 -32.35544   3.648140
#> 5         129       23 -28.75226   3.490686
mine_sequences(tree, engagement, n = 5, which = "expected")
#>   sequence_id n_scored  log_lik perplexity
#> 1         354       25 -6.89971   1.317833
#> 2         365       25 -6.89971   1.317833
#> 3         450       25 -6.89971   1.317833
#> 4         536       25 -6.89971   1.317833
#> 5         539       25 -6.89971   1.317833
```
