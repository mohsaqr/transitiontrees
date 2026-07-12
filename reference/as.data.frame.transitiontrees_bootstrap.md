# Coerce a context tree Bootstrap to a Tidy Data Frame

Uniform tidy-extract: returns the per-pathway summary table
(`object$summary`), so `as.data.frame(boot)` and `summary(boot)` are
interchangeable extractors.

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_bootstrap`.

- row.names, optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame; see
[`bootstrap_pathways`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
for the full column vocabulary.
