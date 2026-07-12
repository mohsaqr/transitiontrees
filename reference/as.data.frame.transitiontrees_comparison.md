# Coerce a context tree Comparison to a Tidy Data Frame

Uniform tidy-extract: returns the per-pathway divergence breakdown
(`object$pathways`) — the consumer-facing detail behind the scalar
`pdist` and the permutation `p_value`.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison`.

- row.names, optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame with columns `pathway`, `count_a`, `count_b`,
`divergence_ab`, `divergence_ba`, `divergence_sym`.
