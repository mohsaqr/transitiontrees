# Implementation notes removed from `?plot_pruning` (2026-11-17)

Still true of R/plot_diagnostics.R; kept for maintainers.

- `.suffix_chain()` builds the chain from the requested context down to the
  root; the requested pathway must be a fitted node (error otherwise; no
  silent fallback to a shorter suffix).
- Decision rule = `prune_tree(criterion = "G2")`: `2N * KL > qchisq(1 -
  alpha, k - 1)`; `pruned` is cumulative from the leaf up (amnesia rule);
  `retained` = below cutoff but a deeper context diverges.
- Opacity: informative 1, retained 0.7, pruned 0.45, root 1. The modal bar
  gets a `grey20` outline. Fill = Okabe-Ito (`.okabe_ito()`), by next state.
- Panels in one row, deepest on the left; y limits `[0, 1.08]`.
- `pathway` check is an unclassed `stopifnot()`; `alpha` and fitted-node
  checks are `transitiontrees_bad_input`.
