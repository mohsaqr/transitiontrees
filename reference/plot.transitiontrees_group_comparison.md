# Plot a Group Comparison

Plot a Group Comparison

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
plot(x, style = c("divergence", "matrix"), top = 15L, alpha = 0.05, ...)
```

## Arguments

- x:

  A `transitiontrees_group_comparison`.

- style:

  One of `"divergence"` (default; top pathways ranked by behavioral JSD,
  significant ones highlighted) or `"matrix"` (heatmap of the
  between-group symmetric-KL distance matrix).

- top:

  Integer. Pathways to show in `"divergence"`. Default 15.

- alpha:

  Numeric. FDR cutoff for highlighting. Default 0.05.

- ...:

  Ignored.

## Value

A ggplot object.
