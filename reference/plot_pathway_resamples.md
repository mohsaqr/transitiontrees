# Plot Bootstrap Resample Distributions per Pathway

Faceted histogram of the bootstrap resample values for a chosen pathway
statistic, one panel per pathway.

## Usage

``` r
plot_pathway_resamples(
  x,
  pathways = NULL,
  stat = c("count", "next_probability", "divergence", "G2"),
  top = 6L,
  bins = 30L
)
```

## Arguments

- x:

  A `transitiontrees_bootstrap` object.

- pathways:

  Character vector of pathway names. `NULL` (default) picks the top
  `top` pathways from the summary.

- stat:

  Character. One of `"count"` (default), `"next_probability"`,
  `"divergence"`, `"G2"`.

- top:

  Integer. Default 6.

- bins:

  Integer. Histogram bins. Default 30.

## Value

A ggplot object.

## Examples

``` r
seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
boot <- bootstrap_pathways(context_tree(seqs, max_depth = 1L),
                           iter = 50L)
plot_pathway_resamples(boot, stat = "count")
```
