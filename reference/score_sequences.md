# Per-Sequence Scoring

Returns one row per held-out sequence with its log-likelihood, number of
scored positions, and per-sequence perplexity
(`exp(-log_lik / n_scored)`).

## Usage

``` r
score_sequences(tree, newdata)
```

## Arguments

- tree:

  A `transitiontrees`.

- newdata:

  Sequence data.

## Value

A data.frame with one row per sequence in `newdata` and columns
`sequence_id` (row / element index of `newdata`), `n_scored` (positions
scored), `log_lik` (natural-log likelihood) and `perplexity`.

## Examples

``` r
set.seed(1)
fit  <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 1L)
new  <- replicate(5, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
score_sequences(tree, new)
#>   sequence_id n_scored   log_lik perplexity
#> 1           1       10 -10.57935   2.880417
#> 2           2       10 -10.54631   2.870915
#> 3           3       10 -11.09127   3.031711
#> 4           4       10 -10.99497   3.002655
#> 5           5       10 -10.71298   2.919165
```
