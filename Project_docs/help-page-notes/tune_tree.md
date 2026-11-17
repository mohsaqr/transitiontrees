# Implementation notes removed from `?tune_tree` and its methods (2026-11-17)

Still true of R/tune.R; kept for maintainers.

- **Folds** (`.pt_make_folds()`): `set.seed(seed)` is called on the global
  RNG and not restored; sequences are shuffled with `sample.int(n)` and
  assigned round-robin, so fold sizes differ by at most one. Every grid
  point reuses the same folds.
- **Input**: long logs go through `.ct_maybe_reshape()` (same path as
  `context_tree()`); then `.ct_traj(.ct_coerce(data))`. The alphabet is
  taken from all sequences before splitting, so every fold's tree shares it.
- **Grid**: `expand.grid(max_depth, min_count, smoothing_i, prune)`.
  Smoothing labels come from `.pt_smoothing_label()`, e.g.
  `floor(ymin=0.001, rule=interpolate)` for the default `"floor"`.
  `.pt_smoothing_grid()` treats a list whose `$method` is set, or whose first
  element is a single known method name, as ONE spec; otherwise as a list of
  specs.
- **Pruning inside CV**: `prune_tree(tr, criterion = "G2", alpha = alpha)`;
  the other criteria are not available in the grid.
- **Scoring**: `logLik.transitiontrees(tr, newdata = test)`; perplexity is
  `exp(-sum(ll) / sum(nobs))` over the folds that scored. A fold that errors
  is caught (`tryCatch`), counted in `folds_failed`, and its message is kept
  for the `transitiontrees_fold_failed` warning (first message only).
- **Selection**: rows ordered by `perplexity` with `na.last = TRUE`; ties keep
  `expand.grid` order. `best` = first row with `folds_failed == 0` and finite
  perplexity; stored as the `"best"` attribute.
- **plot()**: `geom_point(shape = 8, colour = "#D55E00")` for the star; facets
  `prune_label ~ smoothing` when both prune values are present, else
  `~ smoothing`; `theme_minimal(base_size = 11)`.
- The column formerly called `nmin` in the grid is `min_count` (renamed
  2026-11-17).
