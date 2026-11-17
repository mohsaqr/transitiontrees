# Generate sequences from a fitted transition tree

`generate_sequences()` draws new sequences from the process described by
a fitted tree. Each sequence starts from a state drawn from the overall
distribution of states (or from a given first state), and each further
state is drawn from the next-state distribution of the longest context
of the tree that matches the states generated so far. Generated
sequences can be compared with the observed ones to check whether the
tree reproduces features of the data that it was not fitted to directly,
such as run lengths or the time spent in each state.

Draws use the session's random number generator; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first, or use
[`simulate.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/simulate.transitiontrees.md)
with `seed`, for reproducible sequences. It returns one sequence per
row, as a character matrix or, with `tidy = TRUE`, a data frame.

## Usage

``` r
generate_sequences(tree, n = 5L, length = 10L, start = NULL, tidy = FALSE)
```

## Arguments

- tree:

  A `transitiontrees` object.

- n:

  Number of sequences. Default 5.

- length:

  Number of states in each sequence. Default 10.

- start:

  First state of each sequence, a character vector of length `n`. `NULL`
  (default) draws it from the root distribution.

- tidy:

  Return a `data.frame` rather than a matrix. Default `FALSE`.

## Value

With `tidy = TRUE`, a `data.frame` with one row per sequence and columns
`T1` to `T<length>`, which
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
accepts as input. With `tidy = FALSE`, a character matrix with the same
layout.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree, `n` or
`length` is not a single positive number, or `start` does not have
length `n`.

## See also

[`simulate.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/simulate.transitiontrees.md),
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
set.seed(1)
generate_sequences(tree, n = 5, length = 12, tidy = TRUE)
#>           T1         T2      T3         T4         T5         T6         T7
#> 1    Average Disengaged Average    Average     Active     Active     Active
#> 2    Average Disengaged Average    Average    Average    Average    Average
#> 3     Active     Active  Active Disengaged Disengaged Disengaged Disengaged
#> 4 Disengaged    Average Average    Average    Average     Active     Active
#> 5    Average    Average  Active    Average     Active    Average    Average
#>           T8         T9        T10     T11     T12
#> 1     Active     Active     Active  Active  Active
#> 2 Disengaged Disengaged Disengaged Average Average
#> 3 Disengaged    Average    Average Average Average
#> 4    Average    Average     Active  Active Average
#> 5    Average    Average    Average Average Average
```
