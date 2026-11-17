# Implementation notes removed from `?plot_difference` (2026-11-17)

Still true of R/plot_diagnostics.R; kept for maintainers.

- Tile layout: contexts = intersection of the two trees' node names, filtered
  by `depth` and by `n >= min_count` in BOTH trees.
- Residual: `o` = 2 x K matrix of the two nodes' `counts`; `e =
  .cg_expected(o)`; cell = `(o[1, ] - e[1, ]) / sqrt(e[1, ])`, non-finite
  set to 0. Probability measure uses the trees' smoothed `prob`.
- Labels shown when `|value| >= 0.5` (residual) or `>= 0.01` (probability);
  formats `%+.1f` / `%+.2f`. Rows ordered by `sum(abs(value))`.
- Stars: `comparison$pathways$jsd_padj < alpha` appends `" *"` to the label.
- Colours: `scale_fill_gradient2(low = "#2166AC", mid = "white", high =
  "#B2182B")`, symmetric limits.
- Tree layout (`.diff_tree()`): pooled tree refitted with the FIRST group's
  `max_depth` and `nmin` (default smoothing); delta = context count share in
  group 1 minus group 2, each over that group's root count; drawn on the
  horizontal layout; internal nodes labelled with `.pt_last_state()` (oldest
  token), leaves with the full context; subtitle says "phase" (legacy
  wording).
