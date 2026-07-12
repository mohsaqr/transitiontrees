# Log-Likelihood of a context tree

Returns a `logLik` object compatible with
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html),
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html), and the rest of the
model-comparison toolchain. If `newdata` is `NULL`, returns the
in-sample log-likelihood computed from the fitted node counts. Otherwise
returns the held-out log-likelihood scoring `newdata` under the fitted
tree.

## Usage

``` r
# S3 method for class 'transitiontrees'
logLik(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `transitiontrees`.

- newdata:

  Optional. Sequence data in any format accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).
  `NULL` (default) returns in-sample log-likelihood.

- ...:

  Ignored.

## Value

A `logLik` object with attributes `nobs` and `df` (number of free
parameters in the fitted tree).

## Details

Out-of-vocabulary handling: when `newdata` contains a state not in the
tree's alphabet, the transition *into* that state is omitted from
scoring (it is not penalised), and the transition *out* of it is scored
against the root context (the unseen state cannot extend a history). The
reported `nobs` therefore counts only the positions actually scored.
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
and
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
inherit the same behaviour.

## Examples

``` r
# \donttest{
tree <- context_tree(matrix(sample(c("A","B","C"), 200, TRUE), 20),
                     max_depth = 2, min_count = 2)
logLik(tree)
#> 'log Lik.' -210.7309 (df=26)
AIC(tree); BIC(tree)
#> [1] 473.4618
#> [1] 559.2181
# }
```
