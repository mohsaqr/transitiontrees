# Coerce a Group Comparison to a Tidy Data Frame

Uniform tidy-extract: returns the per-pathway comparison table
(`object$pathways`) so `as.data.frame(cmp)` yields the full divergence /
usage breakdown as a base `data.frame`.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_group_comparison`.

- row.names, optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame with one row per tested context; see
[`compare_groups`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for the column vocabulary.

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
as.data.frame(cmp)
#>   pathway depth count_total count_x count_y modal_x modal_y flips  jsd_bits
#> 1       B     1         226     175      51       B       C  TRUE 0.1706917
#> 2 (start)     0         640     320     320       B       C  TRUE 0.1619964
#> 3       A     1         117      49      68       B       C  TRUE 0.1418159
#> 4       C     1         217      56     161       B       C  TRUE 0.1163298
#>   jsd_p jsd_padj   usage_g2 usage_p usage_padj
#> 1  0.01     0.01 119.098496    0.01      0.015
#> 2  0.01     0.01         NA      NA         NA
#> 3  0.01     0.01   3.914325    0.02      0.020
#> 4  0.01     0.01  85.666722    0.01      0.015
```
