# Predict Next-State Probabilities from a Context Tree

Predict Next-State Probabilities from a Context Tree

## Usage

``` r
# S3 method for class 'transitiontrees'
predict(object, newdata, type = c("prob", "class"), ...)
```

## Arguments

- object:

  A `transitiontrees`.

- newdata:

  Either (i) a list of character vectors (each is the "history" leading
  up to the prediction point), (ii) a wide data.frame / matrix whose
  rows are histories, or (iii) a single character vector treated as one
  history.

- type:

  One of `"prob"` (default; named numeric matrix of next-state
  probabilities) or `"class"` (character vector of modal predictions).

- ...:

  Ignored.

## Value

If `type = "prob"`: a matrix with one row per history and one column per
state. A list/data.frame/matrix `newdata` always returns a matrix (1 x k
for a single-history container); a bare character vector returns a named
vector for interactive convenience. If `type = "class"`: a character
vector of modal predictions.

## Examples

``` r
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
predict(tree, newdata = list(c("A","B"), c("C","C","B")))
#>              A         B         C
#> [1,] 0.3333333 0.3518519 0.3148148
#> [2,] 0.4117647 0.2941176 0.2941176
predict(tree, newdata = list(c("A","B")), type = "class")
#> [1] "B"
predict(tree, newdata = c("A","B"))   # bare vector → named vector
#>         A         B         C 
#> 0.3333333 0.3518519 0.3148148 
# }
```
