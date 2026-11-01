# Print a Group Comparison

Print method for a `transitiontrees_group_comparison`: a header with the
groups and permutation setup, the omnibus behavioral / usage statistics,
and the top pathways ranked by behavioral divergence.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
print(x, n = 10L, digits = 3L, ...)
```

## Arguments

- x:

  A `transitiontrees_group_comparison` object.

- n:

  Integer. Number of top pathways to print. Default 10.

- digits:

  Integer. Numeric digits for the printed tables. Default 3.

- ...:

  Ignored.

## Value

`x` invisibly.

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
print(cmp, n = 3L)
#> <transitiontrees_group_comparison>  2 groups, 99 permutations
#>   groups   : x, y
#>   contexts : 4 tested
#> 
#> omnibus (permutation):
#>        axis                 statistic  value p_value
#>  behavioral count-weighted JSD (bits) 184.09    0.01
#>       usage                   sum G^2 208.68    0.01
#> 
#> top pathways by behavioral divergence (JSD):
#>  pathway depth count_total modal_x modal_y flips jsd_bits jsd_padj usage_padj
#>        B     1         226       B       C  TRUE    0.171     0.01      0.015
#>  (start)     0         640       B       C  TRUE    0.162     0.01         NA
#>        A     1         117       B       C  TRUE    0.142     0.01      0.020
#> # ... 1 more (use as.data.frame(x) for the full table)
```
