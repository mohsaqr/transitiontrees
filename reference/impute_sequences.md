# Fill in missing states of incomplete sequences

`impute_sequences()` fills in the missing states of incomplete sequences
from a fitted tree. Sequences recorded over time often have gaps:
positions at which the state was not observed although the sequence
continues afterwards. The tree gives, for each gap, a distribution of
the missing state conditional on the states observed before it, which is
used to fill the gap with the most likely state or with a random draw.

Gaps are filled from left to right. Each gap is predicted from the
longest context of the tree that matches the states before it, and a
filled state becomes part of the history for the gaps that follow. Only
gaps followed by an observed state are filled: missing positions at the
end of a sequence mark where it stops, and a sequence with no observed
state is returned unchanged. The prediction uses past states only; the
states observed after a gap do not inform it.

With `method = "modal"` the result is deterministic and the filled
states are the single most likely ones, which understates the
variability of the process; with `method = "prob"` the states are drawn
from the predicted distribution, and repeated imputations with different
seeds reflect the uncertainty of the filled values. It returns `newdata`
in the same form with its internal gaps filled.

## Usage

``` r
impute_sequences(tree, newdata, method = c("modal", "prob"), seed = NULL)
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Sequences with gaps, marked by `NA` or `""`: a list of character
  vectors, a character matrix or `data.frame` with one sequence per row,
  or a single character vector.

- method:

  `"modal"` (default) fills each gap with the most likely state;
  `"prob"` draws it from the predicted distribution.

- seed:

  Random seed for `method = "prob"`, set locally without changing the
  session's random number stream. `NULL` (default) uses the current
  stream.

## Value

`newdata` in the same form (list, matrix, `data.frame` or vector) with
its internal gaps filled.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree or
`newdata` is not a list, matrix, `data.frame` or character vector.

## See also

[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md),
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
gappy <- list(c("Active", NA, "Active", "Average"),
              c("Disengaged", "Disengaged", NA, NA, "Average", NA))
impute_sequences(tree, gappy)
#> [[1]]
#> [1] "Active"  "Active"  "Active"  "Average"
#> 
#> [[2]]
#> [1] "Disengaged" "Disengaged" "Disengaged" "Disengaged" "Average"   
#> [6] NA          
#> 
impute_sequences(tree, gappy, method = "prob", seed = 1)
#> [[1]]
#> [1] "Active"  "Active"  "Active"  "Average"
#> 
#> [[2]]
#> [1] "Disengaged" "Disengaged" "Disengaged" "Disengaged" "Average"   
#> [6] NA          
#> 
```
