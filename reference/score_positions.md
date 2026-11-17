# Score each position of each sequence under a fitted tree

`score_positions()` returns, for every position of every sequence, the
context of the tree that was used to predict it, the state that
occurred, the probability the tree gave that state, and its
log-likelihood contribution. Low probabilities mark the moves the tree
did not expect, which locates where in a sequence it departs from the
typical pathways.

With `worst`, only the positions with the lowest predicted probability
are returned, in increasing order of probability. Summed by sequence,
the log-likelihood contributions give
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md).

## Usage

``` r
score_positions(tree, newdata, worst = NULL)
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Sequences to score, in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- worst:

  Number of least expected positions to return. `NULL` (default) returns
  all positions, in sequence order.

## Value

A `data.frame` with one row per scored position and columns
`sequence_id`, `position`, `matched_context` (the context used for the
prediction; `"(start)"` for the root), `observed` (the state that
occurred), `predicted_prob` (its probability) and `log_lik` (its natural
logarithm).

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree, `worst` is
not a single positive number, or `newdata` has no usable sequence.

## See also

[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md),
[`plot_predictive()`](https://pak.dynasite.org/transitiontrees/reference/plot_predictive.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
score_positions(tree, engagement, worst = 5)
#>   sequence_id position  matched_context   observed predicted_prob   log_lik
#> 1           1       18 Active -> Active Disengaged     0.02078522 -3.873513
#> 2           2       22 Active -> Active Disengaged     0.02078522 -3.873513
#> 3           3       12 Active -> Active Disengaged     0.02078522 -3.873513
#> 4           4       21 Active -> Active Disengaged     0.02078522 -3.873513
#> 5           5        3 Active -> Active Disengaged     0.02078522 -3.873513
```
