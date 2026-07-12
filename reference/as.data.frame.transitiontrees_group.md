# Coerce a Group of Trees to One Tidy Data Frame

Row-binds each group's
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
table, tagged with a leading `group` column, so the whole batch is one
tidy frame ready to filter, facet, or join.

## Usage

``` r
# S3 method for class 'transitiontrees_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_group`.

- row.names, optional:

  Ignored.

- ...:

  Forwarded to
  [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Value

A data.frame: the canonical pathway columns with a leading `group`
column identifying the source tree.
