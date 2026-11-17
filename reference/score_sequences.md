# Score each sequence under a fitted tree

`score_sequences()` scores each sequence by how well the fitted tree
predicts it: its log-likelihood, the number of positions scored, and its
perplexity (\\\exp(-\ell / n)\\). Sequences with high perplexity follow
pathways the tree finds unlikely, which makes the score useful for
finding atypical sequences or for checking how well a tree fitted on one
group describes another.

Each position is predicted from the longest context of the tree that
matches the preceding states, as in
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md);
states not among the tree's states are not scored. It returns one row
per scored sequence. For the score of each position use
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md);
for the most or least typical sequences directly, use
[`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md).

## Usage

``` r
score_sequences(tree, newdata)
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Sequences to score, in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Value

A `data.frame` with one row per scored sequence and columns
`sequence_id` (position of the sequence among the non-empty sequences of
`newdata`), `n_scored` (positions scored), `log_lik` (natural logarithm)
and `perplexity`.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`newdata` has an unsupported type or no usable sequence.

## See also

[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md),
[`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md),
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
head(score_sequences(tree, engagement))
#>   sequence_id n_scored   log_lik perplexity
#> 1           1       23 -19.32997   2.317371
#> 2           2       23 -17.60598   2.150021
#> 3           3       24 -19.79227   2.281146
#> 4           4       24 -20.34323   2.334119
#> 5           5       24 -21.69346   2.469199
#> 6           6       23 -14.81104   1.904002
```
