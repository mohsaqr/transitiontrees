# Predict the next state after each history

[`predict()`](https://rdrr.io/r/stats/predict.html) predicts the state
that follows each history from the longest context of a fitted tree that
matches the end of that history. Histories longer than the tree's depth
are read only as far back as the tree remembers, and histories whose
recent states never formed a context fall back to a shorter suffix or to
the root. The prediction is therefore always defined, and it uses as
much of the past as the data supported when the tree was fitted.

It returns, for each history, the probability of every next state, or
the most likely next state with `type = "class"`. With `tidy = TRUE` it
returns one row per history with the history, the context that was
matched, the most likely next state and its probability, and the
probability of every state.

To score observed sequences rather than predict from given histories,
use
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md);
to look up a single pathway, use
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md);
to fill in missing states, use
[`impute_sequences()`](https://pak.dynasite.org/transitiontrees/reference/impute_sequences.md).

## Usage

``` r
# S3 method for class 'transitiontrees'
predict(object, newdata, type = c("prob", "class"), tidy = FALSE, ...)
```

## Arguments

- object:

  A `transitiontrees` object.

- newdata:

  Histories, oldest state first: a list of character vectors, a wide
  `data.frame` or matrix with one history per row (missing or empty
  cells are dropped), or a single character vector. A history given as
  one arrow-notation string (`"A -> B"`) is split into its states.

- type:

  `"prob"` (default) for next-state probabilities or `"class"` for the
  most likely next state.

- tidy:

  Return a `data.frame` with one row per history. Default `FALSE`.

- ...:

  Not used.

## Value

With `tidy = TRUE`, a `data.frame` with one row per history and columns
`id` (position in `newdata`), `history` (arrow notation; `"(start)"`
when empty), `context` (the matched context), `predicted` (most likely
next state), `probability` (its probability) and one column per state
with the next-state probabilities. With `tidy = FALSE` and
`type = "prob"`, a numeric matrix with one row per history and one
column per state, or a named vector when `newdata` is a single character
vector; with `type = "class"`, a character vector of most likely next
states.

## Conditions

`transitiontrees_bad_input` when `newdata` is not a list, data frame,
matrix or character vector. A warning of class
`transitiontrees_unknown_state` when a history contains a state outside
the tree's alphabet; the prediction then uses the longest context after
the last unknown state.

## See also

[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md),
[`impute_sequences()`](https://pak.dynasite.org/transitiontrees/reference/impute_sequences.md),
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
histories <- list(c("Active", "Active"),
                  c("Average", "Disengaged", "Disengaged"),
                  character(0))
predict(tree, newdata = histories, tidy = TRUE)
#>   id                             history                             context
#> 1  1                    Active -> Active                    Active -> Active
#> 2  2 Average -> Disengaged -> Disengaged Average -> Disengaged -> Disengaged
#> 3  3                             (start)                             (start)
#>    predicted probability     Active   Average Disengaged
#> 1     Active   0.7852194 0.78521940 0.1939954 0.02078522
#> 2 Disengaged   0.6000000 0.02222222 0.3777778 0.60000000
#> 3    Average   0.4347594 0.37593583 0.4347594 0.18930481
predict(tree, newdata = histories, type = "class")
#> [1] "Active"     "Disengaged" "Average"   
```
