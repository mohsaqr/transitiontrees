# Implementation notes removed from `?impute_sequences`, `?mine_contexts` and `?mine_sequences` (2026-11-17)

Still true of R/impute.R and R/mine.R; kept for maintainers.

- `impute_sequences()`: gaps are `NA` or `""` only (TraMineR `"%"`/`"*"`
  are not treated as gaps here, unlike in `context_tree()`). The history
  for a gap is every non-missing state before it, including states filled
  earlier (sequential loop, justified in the code). Positions after the
  last observed state are left as they are. `seed` is applied with
  `set.seed()` and the caller's `.Random.seed` is restored on exit (or
  removed if it did not exist). A data.frame input is returned as a
  character data.frame with the same dimnames.
- `mine_contexts()`: `is_modal` uses `which.max()` (ties broken by
  alphabet order); `count` is numeric to match the empty template; the
  root is labelled `"(start)"`. `min_prob`/`max_prob` filters are applied
  after `min_count`; sort `order(-prob, -count)`.
- `mine_sequences()` sorts `score_sequences()` output by `-perplexity`
  (surprising) or `perplexity` (expected) and takes `head(n)`.
