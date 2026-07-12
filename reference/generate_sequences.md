# Sample Sequences from a Fitted Context Tree

Sample Sequences from a Fitted Context Tree

## Usage

``` r
generate_sequences(tree, n = 5L, length = 10L, start = NULL)
```

## Arguments

- tree:

  A `transitiontrees`.

- n:

  Integer. Number of sequences to sample.

- length:

  Integer. Length of each sampled sequence.

- start:

  NULL or character vector. If NULL (default), each sequence starts from
  the root marginal; otherwise must have length `n` and supply the first
  state of each sequence.

## Value

A character matrix of dimension `n x length`.

## Examples

``` r
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
generate_sequences(tree, n = 5, length = 10)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] "B"  "B"  "B"  "C"  "B"  "B"  "C"  "C"  "C"  "C"  
#> [2,] "A"  "B"  "C"  "C"  "A"  "A"  "A"  "A"  "C"  "B"  
#> [3,] "A"  "C"  "C"  "A"  "C"  "C"  "A"  "B"  "A"  "C"  
#> [4,] "C"  "C"  "A"  "B"  "B"  "C"  "A"  "C"  "B"  "A"  
#> [5,] "C"  "A"  "A"  "B"  "A"  "B"  "B"  "C"  "A"  "B"  
# }
```
