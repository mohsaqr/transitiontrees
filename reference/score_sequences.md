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

A data.frame with columns `sequence_id`, `n_scored`, `log_lik`,
`perplexity`.

## Examples

``` r
fit  <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 1L)
new  <- replicate(5, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
score_sequences(tree, new)
#>   sequence_id n_scored   log_lik perplexity
#> 1           1       10 -10.34191   2.812829
#> 2           2       10 -11.38681   3.122645
#> 3           3       10 -10.96508   2.993695
#> 4           4       10 -11.80214   3.255072
#> 5           5       10 -11.80038   3.254499
```
