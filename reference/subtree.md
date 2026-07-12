# Extract the Subtree Rooted at a Pathway

Returns a new `transitiontrees` containing only the queried node and its
descendants. Node names are kept absolute (so the original pathway is
preserved as a key), but the returned object has a `local_root`
attribute pointing at the queried pathway.

## Usage

``` r
subtree(tree, pathway)
```

## Arguments

- tree:

  A `transitiontrees`.

- pathway:

  Character. The root pathway (arrow-notation or character vector). Must
  exist in the tree.

## Value

A new `transitiontrees` whose nodes and edges are restricted to
descendants of `pathway`. The alphabet, smoothing, and other
hyperparameters are copied unchanged. Printing the subtree reports the
context it was cut at (also stored in `attr(., "local_root")` for
programmatic use).

## Examples

``` r
# \donttest{
set.seed(1)
m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr <- context_tree(m, max_depth = 2L, min_count = 3L)
subtree(tr, "A")   # prints "subtree of: A" in the banner
#> <transitiontrees>  4 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : A, B, C
#>   fit on   : 20 sequences, 200 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 3
#>   subtree of: A
# }
```
