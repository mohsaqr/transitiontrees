# Perplexity of a context tree

`exp(-mean log-likelihood per observation)`, the standard
language-modelling evaluation metric. Lower is better. A perplexity of
\\k\\ on an alphabet of size \\\|S\|\\ means the model is as predictive
as a uniform distribution over \\k\\ symbols. \\k = \|S\|\\ is the
uniform baseline; \\k = 1\\ is perfect deterministic prediction.

## Usage

``` r
perplexity(tree, newdata = NULL)
```

## Arguments

- tree:

  A `transitiontrees`.

- newdata:

  Sequence data; `NULL` (default) returns in-sample perplexity.

## Value

Numeric scalar.

## Examples

``` r
# \donttest{
tree <- context_tree(matrix(sample(c("A","B","C"), 200, TRUE), 20),
                     max_depth = 2, min_count = 2)
perplexity(tree)
#> [1] 2.715389
# }
```
