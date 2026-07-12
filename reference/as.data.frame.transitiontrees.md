# Coerce a context tree to a Tidy Data Frame

Returns the canonical tidy node table — identical to
`tree_pathways(tree)`. Lets users do `as.data.frame(tree)` and
immediately filter, sort, or export with base-R idioms.

## Usage

``` r
# S3 method for class 'transitiontrees'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees`.

- row.names, optional:

  Ignored.

- ...:

  Forwarded to
  [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).

## Value

A data.frame with columns `pathway`, `depth`, `count`, `likely_next`,
`next_probability`, `divergence`, `changes_prediction`. See
[`tree_pathways`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md).
