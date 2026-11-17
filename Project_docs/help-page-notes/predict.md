# Implementation notes removed from `?predict.transitiontrees`, `?generate_sequences` and `?simulate.transitiontrees` (2026-11-17)

Still true of R/predict.R; kept for maintainers.

- **Context matching** (`.ct_match_context()`): walks k from
  `min(max_depth, length(history))` down to 1 and returns the first
  `paste(tail(history, k), collapse = " -> ")` that is a node; else the
  root. A node missing from the list would yield a uniform distribution
  (only possible if the root itself were missing).
- **predict() output shapes** (tidy = FALSE): a list, data.frame or matrix
  `newdata` always yields a matrix (1 x k for a single history); only a
  bare character vector (no `dim`) collapses to a named vector. A
  character *matrix* is a container. Matrix/data.frame rows drop `NA` and
  `""` cells (all positions, not only trailing). List histories are used
  as given (no NA cleaning).
- `type = "class"` and the tidy `predicted` column break ties with
  `max.col(ties.method = "first")` (alphabet order). With `tidy = TRUE`
  `type` is ignored; state columns go through `make.unique()` in case a
  state name collides with a fixed column name.
- `generate_sequences()`: two nested loops (sequential dependency on the
  generated history); `start` is not checked against the alphabet. Uses
  `sample(alphabet, 1, prob = p)`.
- `simulate.transitiontrees(seed = )` calls `set.seed(seed)` on the global
  RNG and does not restore the previous state (unlike
  `impute_sequences(seed = )`, which restores it).
