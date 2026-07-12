# Mine Sequences by Predictive Surprise

Rank held-out sequences by how well the fitted tree predicts them. A
tidy pattern-mining table: surface the subsequences the model finds most
*surprising* (poor fit, high perplexity) or most *expected* (good fit,
low perplexity).

## Usage

``` r
mine_sequences(tree, newdata, n = 10L, which = c("surprising", "expected"))
```

## Arguments

- tree:

  A `transitiontrees`.

- newdata:

  Sequence data in any format accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- n:

  Integer. Number of sequences to return. Default 10.

- which:

  One of `"surprising"` (default; highest perplexity first) or
  `"expected"` (lowest perplexity first).

## Value

A data.frame with the
[`score_sequences`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md)
columns (`sequence_id`, `n_scored`, `log_lik`, `perplexity`), the top
`n` by the chosen direction.

## Examples

``` r
fit  <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 2L)
new  <- replicate(20, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
mine_sequences(tree, new, n = 5, which = "surprising")
#>   sequence_id n_scored   log_lik perplexity
#> 1           2       10 -12.49314   3.487950
#> 2           7       10 -12.22352   3.395164
#> 3          12       10 -11.91339   3.291486
#> 4           5       10 -11.77453   3.246095
#> 5          11       10 -11.61271   3.193989
```
