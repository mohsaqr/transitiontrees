# Number of Contexts (Nodes) in a Tree

The count of contexts the tree represents — an intuitive accessor for
`length(tree$nodes)` (the number printed in the tree banner).

## Usage

``` r
n_nodes(tree)
```

## Arguments

- tree:

  A `transitiontrees` or `transitiontrees_group`.

## Value

An integer. For a `transitiontrees_group`, a named integer vector with
one count per group.

## Examples

``` r
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
n_nodes(context_tree(seqs, max_depth = 3L))
#> [1] 40
# }
```
