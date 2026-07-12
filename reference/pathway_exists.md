# Test Whether a Pathway Exists in the Tree

Test Whether a Pathway Exists in the Tree

## Usage

``` r
pathway_exists(tree, pathway)
```

## Arguments

- tree:

  A `transitiontrees`.

- pathway:

  Character. Pathway as arrow-notation string or character vector.

## Value

Logical scalar.

## Examples

``` r
# \donttest{
tr <- context_tree(matrix(sample(c("A","B"), 50, TRUE), 5),
                   max_depth = 2L, min_count = 1L)
pathway_exists(tr, "A")
#> [1] TRUE
# }
```
