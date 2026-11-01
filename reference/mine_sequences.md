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
set.seed(1)
fit  <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(fit, max_depth = 2L)
new  <- replicate(20, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
mine_sequences(tree, new, n = 5, which = "surprising")
#>   sequence_id n_scored   log_lik perplexity
#> 1           3       10 -13.12440   3.715229
#> 2          17       10 -12.89721   3.631772
#> 3           1       10 -12.80947   3.600046
#> 4          14       10 -12.25846   3.407048
#> 5          11       10 -11.81596   3.259574
```
