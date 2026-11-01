# Summarise a Group Comparison

Prints a compact verdict for a `transitiontrees_group_comparison`: the
omnibus behavioral and usage permutation p-values, and how many pathways
pass the FDR cutoff on each axis or flip their modal next state between
groups. Returns the per-pathway table invisibly.

## Usage

``` r
# S3 method for class 'transitiontrees_group_comparison'
summary(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `transitiontrees_group_comparison` object.

- alpha:

  Numeric. FDR cutoff used when counting significant pathways. Default
  0.05.

- ...:

  Ignored.

## Value

Invisibly, the per-pathway data.frame (`object$pathways`); see
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
summary(cmp)
#> Group comparison: x vs y
#>   behavioral omnibus p = 0.010 | usage omnibus p = 0.010
#>   pathways with behavioral difference (jsd_padj < 0.05): 4
#>   pathways with usage difference      (usage_padj < 0.05): 3
#>   pathways whose modal next state flips between groups: 4
```
