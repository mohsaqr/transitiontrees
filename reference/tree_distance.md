# Symmetric KL Distance Between Two context trees

Bare-metal scalar: a count-weighted average of per-context symmetric
Kullback-Leibler divergence over the union of pathways present in either
tree. No null distribution.

## Usage

``` r
tree_distance(tree_a, tree_b, symmetric = TRUE)
```

## Arguments

- tree_a, tree_b:

  context trees with matching alphabets.

- symmetric:

  Logical. `TRUE` (default) returns \\0.5(D\_{KL}(A\\B) +
  D\_{KL}(B\\A))\\; `FALSE` returns \\D\_{KL}(A\\B)\\ only.

## Value

Numeric scalar.

## Examples

``` r
# \donttest{
set.seed(1)
m1 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
m2 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr1 <- context_tree(m1, max_depth = 2L, min_count = 3L)
tr2 <- context_tree(m2, max_depth = 2L, min_count = 3L)
tree_distance(tr1, tr2)
#> [1] 0.05750785
# }
```
