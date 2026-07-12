# Summary of a Context Tree

Summary of a Context Tree

## Usage

``` r
# S3 method for class 'transitiontrees'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees`.

- ...:

  Ignored.

## Value

A `summary.transitiontrees` object. The `$table` slot is the canonical
pathway data.frame from
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
(columns `pathway`, `depth`, `count`, `likely_next`, `next_probability`,
`divergence`, `changes_prediction`), re-sorted by `(depth, -count)` so
the structural tree order is read top-to-bottom.
