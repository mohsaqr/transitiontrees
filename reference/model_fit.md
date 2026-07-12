# Model-Fit Scalars in One Call

Bundles the standard goodness-of-fit scalars for a fitted tree into one
tidy row: `logLik`, the parameter count `df`, the observation count
`nobs`, `AIC`, `BIC`, and `perplexity`. A one-call replacement for
[`logLik(); nobs(); AIC(); BIC(); perplexity()`](https://rdrr.io/r/stats/logLik.html).

## Usage

``` r
model_fit(tree, newdata = NULL)
```

## Arguments

- tree:

  A `transitiontrees` or `transitiontrees_group`.

- newdata:

  Optional sequence data. If supplied, the scalars are evaluated on it
  (held-out); if `NULL` (default), in-sample.

## Value

A one-row `data.frame` (one row per group for a `transitiontrees_group`)
with columns `logLik`, `df`, `nobs`, `AIC`, `BIC`, `perplexity`.

## Details

With `newdata`, every scalar is computed **out-of-sample** (`AIC`/`BIC`
use the held-out deviance with the model's training `df`). A
`transitiontrees_group` returns one row per group, tagged with a leading
`group` column.

## See also

[`perplexity`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
[`logLik.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md).

## Examples

``` r
# \donttest{
seqs <- replicate(60, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 2L, min_count = 3L)
model_fit(tree)
#>    logLik df nobs     AIC     BIC perplexity
#> 1 -781.96 26  720 1615.92 1734.98   2.962565
# }
```
