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

  Integer or `NULL`. If given, return only the `worst` positions — those
  with the lowest `predicted_prob` (the moves the model was most
  surprised by). Default `NULL` (all positions, in sequence order).

## Value

A data.frame with columns `sequence_id`, `position`, `matched_context`,
`observed`, `predicted_prob`, `log_lik`.

## Examples

``` r
fit  <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 1L)
new  <- replicate(5, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
score_positions(tree, new, worst = 5L)
#>   sequence_id position matched_context observed predicted_prob   log_lik
#> 1           1        4               B        C       0.245614 -1.403994
#> 2           2       10               B        C       0.245614 -1.403994
#> 3           3        2               B        C       0.245614 -1.403994
#> 4           3        6               B        C       0.245614 -1.403994
#> 5           4        2               B        C       0.245614 -1.403994
```
