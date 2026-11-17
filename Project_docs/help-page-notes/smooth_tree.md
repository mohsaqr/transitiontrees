# Implementation notes removed from `?smooth_tree` and `?compare_smoothing` (2026-11-17)

Still true of R/smoothing.R; kept for maintainers.

- **Order of re-smoothing**: nodes are visited by increasing depth (a
  sequential loop, because each node reads its already re-smoothed
  parent from `new_nodes`); the parent is the nearest ancestor present in
  the node list (`.pt_parent_prob()`), since `min_count` can drop
  intermediate suffixes. The root has no parent; back-off schemes then use
  the uniform distribution.
- **Empty contexts** (`n == 0`): floor returns the parent's distribution
  (uniform at the root); the back-off schemes return the parent.
- **floor, `ymin <= 0`**: raw proportions are returned.
- **Kneser-Ney**: the canonical continuation-count formulation needs
  per-state type counts, which the package does not track; the parent's
  smoothed probability is used as the back-off distribution. The old help
  page attributed this approximation to Begleiter, El-Yaniv & Yona (2004),
  Section 3; that section-level attribution was not re-verified and was
  dropped from the help page. The back-off weight is computed as the
  residual `1 - sum(high)` so the result sums to exactly 1, which also
  handles fractional (weighted) counts below `discount`.
- **Hyperparameter validation** (`.pt_validate_smoothing()`): ymin finite
  >= 0; rule in interpolate/cap; alpha finite >= 0; discount and lambda in
  [0, 1].
- **Smoothing spec shapes**: a string, or a list whose `method` entry or
  first element names the method; other named entries override the
  defaults. Both shapes are part of the API contract (CLAUDE.md).
- **Group input**: `smooth_tree()` on a `transitiontrees_group` re-smooths
  each member and keeps the class and `group` attribute.
- **PST parity**: only `floor` with `rule = "interpolate"` is at parity
  with PST; the other four schemes and `rule = "cap"` are the package's
  own designs.
- **compare_smoothing()** calls `n_nodes()` and `perplexity()` on each fit;
  the perplexity is in-sample.
- **Reference correction (2026-11-17)**: Chen & Goodman (1999) pages are
  359-393 (Crossref), not 359-394.
