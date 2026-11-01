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

## Examples

``` r
set.seed(1)
gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                           prob = c(.2,.6,.2)), simplify = FALSE)
gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                           prob = c(.2,.2,.6)), simplify = FALSE)
grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                    max_depth = 1L)
cmp <- compare_groups(grp, iter = 99L)
plot(cmp)

plot(cmp, style = "matrix")
```
