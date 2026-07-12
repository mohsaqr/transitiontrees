# Per-Context Next-State Distributions

Small-multiples bar chart of the full next-state distribution
`P(next | context)`, one panel per context. A per-context probability
display: where
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md)
renders the same numbers as a heatmap, this shows each context's
distribution as its own panel, with the modal bar highlighted.

## Usage

``` r
plot_distributions(tree, contexts = NULL, top = 12L, min_count = 1L)
```

## Arguments

- tree:

  A `transitiontrees`.

- contexts:

  Character vector of pathway strings to show, or `NULL` (default) to
  take the `top` most frequent contexts.

- top:

  Integer. Number of contexts to show when `contexts` is `NULL`. Default
  12.

- min_count:

  Integer. Drop contexts below this count. Default 1.

## Value

A ggplot object.

## Examples

``` r
# \donttest{
seqs <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 2L)
plot_distributions(tree, top = 9)

# }
```
