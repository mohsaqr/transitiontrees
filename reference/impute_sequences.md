# Impute Missing States in Sequences

Fill the gaps in incomplete sequences using a fitted context tree: each
missing state is predicted from the longest matching context of the
states that precede it. Filling proceeds left to right, so a
just-imputed state becomes part of the context for later gaps.

## Usage

``` r
impute_sequences(tree, newdata, method = c("modal", "prob"), seed = NULL)
```

## Arguments

- tree:

  A `transitiontrees`.

- newdata:

  Sequences with gaps: a list of character vectors, a character matrix /
  data.frame (one row per sequence, `NA` or `""` marking a gap), or a
  single character vector.

- method:

  One of `"modal"` (default; fill with the most likely state) or
  `"prob"` (sample from the predicted distribution).

- seed:

  Integer or `NULL`. Optional RNG seed, used only when
  `method = "prob"`.

## Value

The same container shape as `newdata` (list, matrix, data.frame, or
character vector) with internal gaps filled.

## Details

Only **internal** gaps are imputed. A run of trailing `NA` / `""` cells
(end-of-sequence padding in a wide frame) is left untouched, since there
is no observed state after it to mark the sequence as continuing. A
sequence that is entirely missing is returned unchanged (there is
nothing to condition on).

## Examples

``` r
seqs <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 2L)
gappy <- list(c("A", NA, "C"), c("B", "B", NA, "A"))
impute_sequences(tree, gappy)
#> [[1]]
#> [1] "A" "C" "C"
#> 
#> [[2]]
#> [1] "B" "B" "C" "A"
#> 
```
