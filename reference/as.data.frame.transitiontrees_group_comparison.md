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

A data.frame of per-pathway results; see
[`compare_groups`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for the column vocabulary.
