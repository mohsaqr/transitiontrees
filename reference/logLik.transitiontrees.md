# Compute the log-likelihood of a transition tree

[`logLik()`](https://rdrr.io/r/stats/logLik.html) computes the
log-likelihood of a fitted tree: the sum, over every observed position,
of the log probability that the tree assigned to the state that
occurred, given the longest context of the tree that matches the
preceding history. It is the basis of the information criteria and of
perplexity, and allows trees with different depths, pruning or smoothing
to be compared on the same data.

Without `newdata`, the log-likelihood is computed in-sample from the
fitted counts. With `newdata`, the sequences are scored under the fitted
tree, which gives the held-out log-likelihood. A state that is not among
the tree's states is not scored, and the position after it is predicted
from the root; the reported number of observations counts only scored
positions. The number of parameters is the number of contexts times
\\k - 1\\, for \\k\\ states.

The result is a `logLik` object and works with
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html).
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
returns log-likelihood, AIC, BIC and perplexity together as a table.

## Usage

``` r
# S3 method for class 'transitiontrees'
logLik(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `transitiontrees` object.

- newdata:

  Sequences to score, in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).
  `NULL` (default) gives the in-sample log-likelihood.

- ...:

  Not used.

## Value

An object of class `logLik`: the log-likelihood (natural logarithm),
with attributes `df` (number of parameters) and `nobs` (number of scored
positions).

## Conditions

`transitiontrees_bad_input` when `newdata` has an unsupported type or no
usable sequence.

## See also

[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md),
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
[`nobs.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/nobs.transitiontrees.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
logLik(tree)
#> 'log Lik.' -1565.772 (df=26)
AIC(tree)
#> [1] 3183.543
logLik(tree, newdata = engagement)
#> 'log Lik.' -18657.84 (df=26)
```
