# Re-Smooth a Fitted context tree

Replaces every node's probability vector with a new smoothing scheme
without refitting the tree. Walks nodes top-down by depth so each node's
parent is re-smoothed before its children read it.

## Usage

``` r
smooth_tree(tree, smoothing = "floor")
```

## Arguments

- tree:

  A `transitiontrees`.

- smoothing:

  Smoothing specification: either a method name as a string (uses
  defaults for that method's hyperparameters) or a list of the form
  `list(method, ...kwargs)` for explicit hyperparameters. Available
  methods: `"floor"` (`ymin = 0.001`), `"laplace"` (`alpha = 1`),
  `"kneser_ney"` (`discount = 0.75`), `"witten_bell"`,
  `"jelinek_mercer"` (`lambda = 0.5`).

## Value

A new `transitiontrees` with re-smoothed probabilities. Counts and
topology are unchanged.

## Details

For `"kneser_ney"` the canonical continuation-distribution formulation
requires per-state *type counts*. transitiontrees does not track these;
the implementation uses the parent's smoothed probability as the
back-off distribution, an approximation discussed in Begleiter, El-Yaniv
& Yona (2004), *JAIR* 22, §3.

## Examples

``` r
# \donttest{
set.seed(1)
m  <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr <- context_tree(m, max_depth = 2L, min_count = 3L)
smooth_tree(tr, "kneser_ney")
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : A, B, C
#>   fit on   : 20 sequences, 200 observations
#>   smoothing: kneser_ney(discount=0.75)   min_count = 3
#> (start)   n=200    -> B (0.34)
#> |-- A         n=58     -> C (0.43)
#> |   |-- A         n=13     -> C (0.56)
#> |   |-- B         n=24     -> C (0.47)
#> |   `-- C         n=14     -> B (0.50)
#> |-- B         n=63     -> A (0.43)
#> |   |-- A         n=17     -> A (0.54)
#> |   |-- B         n=17     -> B (0.41)
#> |   `-- C         n=22     -> A (0.51)
#> `-- C         n=59     -> B (0.42)
#>     |-- A         n=22     -> B (0.42)
#>     |-- B         n=16     -> C (0.49)
#>     `-- C         n=15     -> B (0.48) 
smooth_tree(tr, list("kneser_ney", discount = 0.5))
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : A, B, C
#>   fit on   : 20 sequences, 200 observations
#>   smoothing: kneser_ney(discount=0.5)   min_count = 3
#> (start)   n=200    -> B (0.34)
#> |-- A         n=58     -> C (0.43)
#> |   |-- A         n=13     -> C (0.55)
#> |   |-- B         n=24     -> C (0.46)
#> |   `-- C         n=14     -> B (0.50)
#> |-- B         n=63     -> A (0.43)
#> |   |-- A         n=17     -> A (0.54)
#> |   |-- B         n=17     -> B (0.41)
#> |   `-- C         n=22     -> A (0.51)
#> `-- C         n=59     -> B (0.42)
#>     |-- A         n=22     -> B (0.42)
#>     |-- B         n=16     -> C (0.50)
#>     `-- C         n=15     -> B (0.48) 
# }
```
