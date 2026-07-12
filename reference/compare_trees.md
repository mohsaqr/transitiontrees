# Compare Two context trees by Symmetric Divergence with Permutation Test

Computes the count-weighted symmetric Kullback-Leibler divergence
between two fitted transitiontrees, then provides a reference
distribution by permuting sequence-to-tree assignments.

Use this to ask: do two cohorts (group A vs. group B, baseline vs.
intervention) generate significantly different pathway dynamics?

## Usage

``` r
compare_trees(tree_a, tree_b = NULL, iter = 200L, seed = 1L, symmetric = TRUE)
```

## Arguments

- tree_a, tree_b:

  context trees fit on data subsets A and B. Alternatively, pass a
  two-element `transitiontrees_group` (from
  `context_tree(..., group =)`) as `tree_a` and leave `tree_b = NULL`;
  its two trees are compared in key order.

- iter:

  Integer. Number of permutations. Default 200.

- seed:

  Integer. RNG seed for reproducibility. Default 1.

- symmetric:

  Logical. Default `TRUE`.

## Value

A `transitiontrees_comparison` object with components:

- pdist:

  observed scalar distance

- null_dist:

  numeric vector, length `iter`

- p_value:

  one-sided p-value (proportion of null at least as extreme as observed)

- pathways:

  per-pathway breakdown data.frame

## Examples

``` r
# \donttest{
set.seed(1)
m1 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
m2 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr1 <- context_tree(m1, max_depth = 2L, min_count = 3L)
tr2 <- context_tree(m2, max_depth = 2L, min_count = 3L)
compare_trees(tr1, tr2, iter = 50)
#> <transitiontrees_comparison>  iter = 50
#>   observed distance : 0.0575
#>   null mean         : 0.0533
#>   p-value           : 0.314
#> 
#> top divergent pathways:
#>  pathway count_a count_b divergence_ab divergence_ba divergence_sym
#>   B -> B      17      12         0.465         0.372          0.418
#>   C -> A      14       9         0.273         0.293          0.283
#>   A -> A      13      29         0.207         0.250          0.228
#>   C -> C      15      21         0.162         0.199          0.181
#>        A      58      72         0.074         0.080          0.077
#>   A -> B      17      17         0.067         0.066          0.066
# }
```
