# Simulate sequences from a fitted transition tree

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws `nsim`
sequences from a fitted tree as
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
does, and sets the random seed first when `seed` is given. It is the
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for fitted
trees.

## Usage

``` r
# S3 method for class 'transitiontrees'
simulate(
  object,
  nsim = 5L,
  seed = NULL,
  length = 10L,
  start = NULL,
  tidy = FALSE,
  ...
)
```

## Arguments

- object:

  A `transitiontrees` object.

- nsim:

  Number of sequences. Default 5.

- seed:

  Random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html), or `NULL`
  (default) to use the current state of the generator.

- length:

  Number of states in each sequence. Default 10.

- start:

  First state of each sequence; see
  [`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).
  Default `NULL`.

- tidy:

  Return a `data.frame` rather than a matrix. Default `FALSE`.

- ...:

  Not used.

## Value

As
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md):
with `tidy = TRUE`, a `data.frame` with one row per sequence and columns
`T1` to `T<length>`; otherwise a character matrix with the same layout.

## Conditions

As
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).

## See also

[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
simulate(tree, nsim = 3, seed = 42, length = 8, tidy = TRUE)
#>           T1      T2         T3         T4      T5         T6         T7
#> 1 Disengaged Average Disengaged Disengaged  Active Disengaged Disengaged
#> 2 Disengaged Average    Average    Average Average Disengaged Disengaged
#> 3    Average Average    Average     Active  Active     Active    Average
#>           T8
#> 1 Disengaged
#> 2     Active
#> 3 Disengaged
```
