# Implementation notes removed from `?plot_predictive` and `?plot_distributions` (2026-11-17)

Still true of R/plot_diagnostics.R; kept for maintainers.

- plot_predictive(): all three types are built from `score_positions(tree,
  newdata)`. "position": grey lines per sequence, points coloured
  `scale_colour_gradient("#D55E00", "#0072B2")`. "logloss": `-log_lik /
  log(2)`; points jittered (width 0.18); per-position mean via
  `aggregate(logloss ~ position)`; dashed overall mean in `#D55E00`.
  "ecdf": `stat_ecdf(geom = "step")`.
- plot_distributions(): contexts from `tree_pathways(min_count, sort_by =
  "count")`, `head(top)` or matched to `contexts`; probabilities are the
  smoothed node `prob`; modal bar `#08519c`, others `#9ecae1`; y limits
  [0, 1].
