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
