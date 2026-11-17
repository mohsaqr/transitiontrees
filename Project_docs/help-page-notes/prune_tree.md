# Implementation notes removed from `?prune_tree` and `?compare_pruning` (2026-11-17)

Still true of R/prune.R; kept for maintainers.

- `prune_tree()` was renamed from `prune()` to avoid collision with other
  `prune()` generics.
- **Algorithm**: repeat { leaves = nodes that are nobody's parent, minus
  the root; for each leaf whose parent is present, evaluate the criterion;
  drop the failing leaves together } until no leaf is dropped. Leaves whose
  parent is absent from the node list are skipped (kept).
- **G2** (`.ct_g2()`): `2 * sum(N * log(N / (n * p_parent)))` over positive
  child counts; expected counts floored at `.Machine$double.eps`; critical
  value `qchisq(1 - alpha, df = k - 1)`. The parent distribution is its
  *smoothed* probability vector.
- **KL** (`.ct_kl()`): in nats; returns `Inf` when the child puts mass on
  a state the parent gives probability 0 (such a leaf is always kept).
- **AIC/BIC**: child log-likelihood uses the child's own smoothed `prob`,
  parent log-likelihood scores the child's counts under the parent's
  `prob`; probabilities floored at `.Machine$double.eps` before `log`.
  AIC penalty `2(k-1)`, BIC penalty `log(n_leaf)(k-1)`.
- `alpha` is used only by G2 and `threshold` only by KL, but both are
  validated for every criterion.
- **Result**: `pruned = TRUE`, `pruning = list(criterion, alpha,
  threshold)`; edges to dropped children removed. Group input keeps the
  class and `group` attribute.
- `compare_pruning()`: `n_nodes` counts the root; `reduction_pct` is
  rounded to one decimal.
- Pruning is not at parity with PST's `prune()` (package's own design;
  see PARITY.md).
- **Reference fix (2026-11-17)**: Ron et al. (1996) citation completed
  (full title, issue 2-3, DOI); DOIs added to Akaike (1974) and Schwarz
  (1978), both verified on Crossref.
