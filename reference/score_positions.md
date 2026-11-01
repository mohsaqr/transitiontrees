# Per-Position Scoring

Returns one row per held-out (sequence, position) with the matched
context, predicted probability of the observed next state, and
log-likelihood contribution. Useful for diagnostic plots showing where
the model is confident vs. surprised.

## Usage

``` r
score_positions(tree, newdata, worst = NULL)
```

## Arguments

- tree:

  A `transitiontrees`.

- newdata:

  Sequence data.

- worst:

  Integer or `NULL`. If given, return only the `worst` positions – those
  with the lowest `predicted_prob` (the moves the model was most
  surprised by), sorted ascending by `predicted_prob`. Must be a single
  positive integer. Default `NULL` (all positions, in sequence order).

## Value

A data.frame with one row per scored (sequence, position) and columns
`sequence_id`, `position`, `matched_context` (the longest context in the
tree that matched the history), `observed` (the state that occurred),
`predicted_prob` (its predicted probability) and `log_lik` (its
natural-log contribution).

## Examples

``` r
set.seed(1)
fit  <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 1L)
new  <- replicate(5, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
score_positions(tree, new, worst = 5L)
#>   sequence_id position matched_context observed predicted_prob   log_lik
#> 1           3        2               B        C      0.2695652 -1.310945
#> 2           4        5               B        C      0.2695652 -1.310945
#> 3           5        9               B        C      0.2695652 -1.310945
#> 4           1        9               A        B      0.2857143 -1.252763
#> 5           2        2               A        B      0.2857143 -1.252763
```
