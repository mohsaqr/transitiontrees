# Simulate Sequences from a Fitted context tree

S3 [`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
`transitiontrees` objects. Wraps
[`generate_sequences`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
with the standard `nsim` argument name and an optional `seed` (set via
[`set.seed()`](https://rdrr.io/r/base/Random.html) when supplied).

## Usage

``` r
# S3 method for class 'transitiontrees'
simulate(object, nsim = 5L, seed = NULL, length = 10L, start = NULL, ...)
```

## Arguments

- object:

  A `transitiontrees`.

- nsim:

  Integer. Number of sequences to simulate. Default 5.

- seed:

  Integer or `NULL`. Optional RNG seed.

- length:

  Integer. Length of each simulated sequence.

- start:

  NULL or character vector. Optional first state for each sequence; see
  [`generate_sequences`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).

- ...:

  Ignored.

## Value

A character matrix of dimension `nsim x length`.
