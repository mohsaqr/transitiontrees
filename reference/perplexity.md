# Compute the perplexity of a transition tree

`perplexity()` computes the perplexity of a fitted tree: the exponential
of the negative mean log-likelihood per observed position, \\\exp(-\ell
/ n)\\. It reads as an effective number of equally likely next states: a
perplexity of \\k\\ for \\k\\ states means the tree predicts no better
than chance, and a perplexity of 1 means every next state is predicted
with certainty. Lower values indicate better prediction.

In-sample perplexity (without `newdata`) favors deep, lightly smoothed
trees that reproduce the fitting data; held-out perplexity, computed on
sequences not used for fitting, measures how well the tree generalizes.
It returns a single number.
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
uses cross-validated perplexity to choose a tree's settings.

## Usage

``` r
perplexity(tree, newdata = NULL)
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Sequences to score, in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).
  `NULL` (default) gives in-sample perplexity.

## Value

A single number; `NA` when no position could be scored.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`newdata` has no usable sequence.

## See also

[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md),
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`plot_predictive()`](https://pak.dynasite.org/transitiontrees/reference/plot_predictive.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
perplexity(tree)
#> [1] 2.310147
perplexity(tree, newdata = engagement)
#> [1] 2.137931
```
