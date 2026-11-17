# Implementation notes removed from `?logLik.transitiontrees`, `?nobs.transitiontrees`, `?perplexity`, `?score_sequences`, `?score_positions`, `?model_fit` and `?n_nodes` (2026-11-17)

Still true of R/likelihood.R; kept for maintainers.

- **In-sample log-likelihood** (`.pt_loglik_in_sample()`): each token is
  attributed to exactly one node, its deepest matching context, via
  "unique counts" (node counts minus the summed counts of its children);
  the sum of unique counts equals `n_obs`. Contributions with zero count
  or zero probability are skipped. Matches PST's `logLik` to 1e-10
  relative (PARITY.md).
- **Held-out scoring** (`.pt_score_walk()`): newdata is coerced with
  `.ct_coerce()`/`.ct_traj()` (same cleaning as the fit); each position's
  context is the longest match of the preceding states, capped at the
  tree's `max_depth`. Out-of-vocabulary observed states are dropped from
  scoring (not penalised); the transition *out of* an unseen state is
  scored against the root (the unseen state cannot extend a history). A
  zero predicted probability gives `log_lik = -Inf`.
- `df = length(nodes) * (k - 1)`; `nobs` attribute = scored positions.
- `nobs()` returns `as.integer(n_obs)` (weighted token total, truncated to
  integer for weighted fits).
- `perplexity()` returns `NA` when `nobs == 0`.
- `score_sequences()`: `sequence_id` is the index among the trajectories
  that survive cleaning (`.ct_traj()`), which equals the row index when
  no row is entirely empty; sequences with no scored position are absent.
- `score_positions(worst = )` sorts by `predicted_prob` ascending (stable).
- `model_fit()`: with zero scored positions every statistic is `NA`
  (`df` kept, `nobs = 0`); AIC = `-2 ll + 2 df`, BIC = `-2 ll + log(n) df`.
  With `newdata`, AIC/BIC use the held-out log-likelihood and the training
  `df`.
- `n_nodes()` counts the root.
