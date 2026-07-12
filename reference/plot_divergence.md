# Lollipop Chart of Pathway Divergence

Lollipop chart of per-pathway Kullback-Leibler divergence from the
(k-1)-suffix. Point size is proportional to pathway count; orange points
mark pathways whose modal next state flips between orders. Annotates
each flip with the prediction change, e.g. "Disengaged -\> Active".

## Usage

``` r
plot_divergence(tree, top = 15L, min_count = 5L, title = NULL, ...)
```

## Arguments

- tree:

  A `transitiontrees`.

- top:

  Integer. Number of pathways to show. Default 15.

- min_count:

  Integer. Drop pathways below this count. Default 5.

- title:

  Character. Plot title; if `NULL` (default) a default title is used.

- ...:

  Ignored.

## Value

A ggplot object.

## Examples

``` r
# \donttest{
seqs <- replicate(50, sample(c("A", "B", "C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3L)
plot_divergence(tree)

# }
```
