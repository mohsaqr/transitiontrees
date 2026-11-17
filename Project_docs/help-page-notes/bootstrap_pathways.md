# Implementation notes removed from `?bootstrap_pathways` and its methods (2026-11-17)

Still true of R/bootstrap.R; kept for maintainers.

- Argument checks (`iter >= 2`, `ci_level`, `alpha`, thresholds in (0, 1),
  `consistency_range` two increasing positives) are named `stopifnot()`
  calls: unclassed errors.
- Counts per depth are precomputed once (`.pt_precompute_for_boot()`); each
  resample sums per-sequence count matrices for `sample.int(n, n, replace =
  TRUE)` (`.pt_resample_count_matrices()`). The loop over `iter` is a
  `for` loop with optional `txtProgressBar`.
- Statistics per pathway (`.pt_pathway_stats_from_counts()`): raw MLE, no
  smoothing; `next_probability` = max count / total; `divergence` =
  KL(child || parent) in bits, `Inf` if a child state has zero parent
  probability; `G2` = `.ct_g2()` (natural log); `changes_prediction` = modal
  index differs from parent's. All `NA` when the pathway is absent in the
  resample (count 0), so `mean_*`/`sd_*`/CIs of these use `na.rm = TRUE`.
- Stability band: `[min, max]` of `orig * consistency_range`; a resample
  with `NA` statistic is outside the band. `p_stability = (n_out + 1) /
  (iter + 1)`; `stable = p_stability < 1 - stability_threshold`.
- `informative_rate = n(G2 > qchisq(1 - alpha, k - 1)) / iter` (non-finite
  G2 not counted, but the denominator is `iter`); `NA` for the root.
- Summary ordering: trust score (`stable + informative`), then
  `stability_rate`, then `count`, all descending.
- Object fields: `summary`, `pathways_orig`, `M_count`,
  `M_next_probability`, `M_divergence`, `M_G2`, `M_changes_prediction`
  (each `iter x n_pathways`, NULL when `keep_resamples = FALSE`), `iter`,
  `stat`, `consistency_range`, `stability_threshold`,
  `informative_threshold`, `alpha_g2`, `ci_level`, `g2_critical_value`,
  `seed`.
- `plot_pathway_resamples()` raises `transitiontrees_bad_input` (not
  `transitiontrees_no_data`) when the resamples were not kept, unlike
  `as.data.frame(x, what = "resamples")`.
- plot(): colours `#0072B2` (stable + informative), `#888888`, `#D55E00`,
  `#cccccc`; `geom_errorbar(height = 0.25)` triggers the ggplot2 4.0 message
  "`height` was translated to `width`".
- Reference history: the earlier text cited "Saqr, Tikka & López-Pernas";
  corrected to the six-author LAK '25 paper (doi:10.1145/3706468.3706513).
