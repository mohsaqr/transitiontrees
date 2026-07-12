# Query the Probability of a Specific Pathway -\> Next State

Returns the probability the fitted tree assigns to a given pathway /
next-state pair. Two lookup modes:

- `exact = TRUE`: the pathway must appear as a node; otherwise returns
  `NA`.

- `exact = FALSE` (default): if the pathway is missing, falls back to
  the longest matching suffix that \*is\* in the tree (mirrors
  [`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md)).

## Usage

``` r
query_pathway(tree, pathway, next_state = NULL, exact = FALSE)
```

## Arguments

- tree:

  A `transitiontrees`.

- pathway:

  Character. The conditioning pathway, either as a single arrow-notation
  string ("A -\> B -\> C") or as a character vector of states
  (`c("A","B","C")`).

- next_state:

  Character. Next-state symbol to query, or `NULL` (default) to return
  the full conditional distribution.

- exact:

  Logical. Default `FALSE` — fall back to longest matching suffix.

## Value

If `next_state` is supplied, a numeric scalar. Otherwise a named numeric
vector indexed by alphabet.

## Examples

``` r
# \donttest{
set.seed(1)
m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr <- context_tree(m, max_depth = 2L, min_count = 3L)
query_pathway(tr, c("A","B"))
#>         A         B         C 
#> 0.5294118 0.2352941 0.2352941 
query_pathway(tr, "A -> B", next_state = "C")
#> [1] 0.2352941
# }
```
