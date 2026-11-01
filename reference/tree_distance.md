# Symmetric KL Distance Between Two context trees

A single number: the count-weighted average of per-context symmetric
Kullback-Leibler divergence (in nats) over the union of pathways present
in either tree. A context missing from one tree is read from its longest
matching suffix there. No null distribution; use
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for a permutation test.

## Usage

``` r
tree_distance(tree_a, tree_b, symmetric = TRUE)
```

## Arguments

- tree_a, tree_b:

  `transitiontrees` objects with matching alphabets.

- symmetric:

  Logical. `TRUE` (default) returns \\0.5(D\_{KL}(A\\B) +
  D\_{KL}(B\\A))\\; `FALSE` returns \\D\_{KL}(A\\B)\\ only.

## Value

Numeric scalar (nats); 0 for identical predictions.

## Examples

``` r
set.seed(1)
m1 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
m2 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr1 <- context_tree(m1, max_depth = 2L, min_count = 3L)
tr2 <- context_tree(m2, max_depth = 2L, min_count = 3L)
tree_distance(tr1, tr2)
#> [1] 0.05750785
```
