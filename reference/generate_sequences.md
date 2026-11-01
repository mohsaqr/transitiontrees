# Sample Sequences from a Fitted Context Tree

Draws sequences by walking the fitted conditional distributions: each
next state is sampled from the longest context in the tree that matches
the history so far. Uses the session RNG, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first (or use
[`simulate.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/simulate.transitiontrees.md)
with `seed`) for reproducible draws.

## Usage

``` r
generate_sequences(tree, n = 5L, length = 10L, start = NULL)
```

## Arguments

- tree:

  A `transitiontrees`.

- n:

  Integer. Number of sequences to sample. Default 5.

- length:

  Integer. Length of each sampled sequence. Default 10.

- start:

  NULL or character vector. If NULL (default), each sequence starts from
  the root marginal; otherwise must have length `n` and supply the first
  state of each sequence.

## Value

A character matrix with `n` rows (one sampled sequence each) and
`length` columns (time-steps).

## Examples

``` r
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
generate_sequences(tree, n = 5, length = 10)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] "A"  "A"  "A"  "A"  "C"  "C"  "B"  "C"  "C"  "A"  
#> [2,] "B"  "C"  "B"  "A"  "B"  "A"  "A"  "B"  "C"  "C"  
#> [3,] "C"  "A"  "C"  "C"  "B"  "C"  "C"  "B"  "C"  "C"  
#> [4,] "A"  "A"  "A"  "A"  "B"  "C"  "B"  "A"  "C"  "B"  
#> [5,] "C"  "B"  "A"  "C"  "C"  "B"  "B"  "C"  "B"  "A"  
```
