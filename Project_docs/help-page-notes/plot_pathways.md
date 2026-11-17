# Implementation notes removed from `?plot_pathways` and `?plot_divergence` (2026-11-17)

Still true of R/plot_pathways.R; kept for maintainers.

- The old `?plot_pathways` said "a side strip on the left encodes the
  pathway count on a log scale". No such strip is drawn; the sentence was
  removed.
- plot_pathways(): rows from `tree_pathways(min_count, sort_by,
  decreasing = TRUE)`, root removed, `head(top)`. Modal cell text white
  bold size 3.2; other cells labelled when prob >= 0.05. Fill
  `scale_fill_gradient("#deebf7", "#08519c")`, limits [0, 1]. Row labels
  monospace; flips prefixed `"> "`. Subtitle prints `n_min = <min_count>`.
- plot_divergence(): rows from `divergent_pathways(top, min_count)`.
  Colour `#0072B2` (same modal) / `#D55E00` (flips); flips annotated
  `"<pct>% <likely_next>"` (the old help page gave the format as
  "Disengaged -> Active", which is not what is drawn). x expansion
  `c(0.02, 0.30)`.
