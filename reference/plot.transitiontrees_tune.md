# Plot a context tree CV Grid

Visualises the held-out perplexity surface returned by
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).
Lines track perplexity vs. `max_depth`; facets split by smoothing scheme
and `prune`; colour encodes `nmin`. The minimum-perplexity configuration
is highlighted with a star.

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
plot(x, ...)
```

## Arguments

- x:

  A `transitiontrees_tune` object.

- ...:

  Ignored.

## Value

A ggplot object.

## Examples

``` r
set.seed(1)
m  <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
tg <- tune_tree(m, max_depth = 1:2, min_count = 5L, prune = FALSE,
                folds = 3L)
plot(tg)
```
