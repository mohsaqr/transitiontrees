# Implementation notes removed from `?plot.transitiontrees` (2026-11-17)

Still true of R/plot.R, R/plot_icicle.R and R/plot_interactive.R; kept for
maintainers.

- **Palette**: `.pt_state_palette()` uses `grDevices::hcl.colors()` "Set 2"
  (<= 8 states), "Set 3" (<= 12) or "Dark 3". `plot_pruning()` uses
  `.okabe_ito()` instead.
- **Fill state**: horizontal = rightmost token (most recent move);
  dendrogram, icicle and interactive = `.pt_last_state()`, which returns the
  FIRST (leftmost, oldest) token despite its name, i.e. the token the node's
  depth adds.
- **horizontal**: leaves ranked top-to-bottom by DFS, evenly spaced on
  [0, 1]; internal y = mean of children's y. Edges are 40-point cosine
  smoothstep curves (`.ct_horizontal_edge_paths()`). Label y offset =
  `0.012 + size_aes * 0.0035` with `size_aes` the sqrt-rescaled count, so
  large hubs push their labels further down. Labels `family = "mono"`,
  size 2.5. Root: black point, size `max(point_size_range) * 0.85`, labelled
  `.ROOT_LABEL`. x expansion `c(0.17, 0.32)`.
- **dendrogram**: leaves evenly spaced in angle `2 * pi * (i - 0.5) / n`;
  internal angle = mean of children's; r = depth. Edges are 24-point arcs
  at the parent radius then a radial segment. Degenerate trees (<= 1 leaf)
  put everything at angle 0. Root drawn last, size
  `max(point_size_range) * 0.8`.
- **icicle**: arc weight = sum of descendant LEAF counts (not the node's own
  count, which can exceed it when sequences end early). Rings:
  thickness 0.68, gap 0.10, ring d spans `d * 0.78 + [0, 0.68]`. Leaf outer
  radius jittered by `runif(-jitter, jitter)` with `jitter_seed`; global RNG
  saved and restored. Polygons sampled at ~1 point per degree. White radial
  separators of width `gap_size` at depth-1 boundaries. Label policy:
  `arc_fraction < label_min_fraction` none; `< label_abbrev_fraction`
  abbreviation (`abbrev_fun` or `substr(, 1, abbrev_len)`, computed over the
  whole alphabet); otherwise the full state name. The in-code comment saying
  `label_abbrev_fraction` defaults to NULL is stale (default 0.08; NULL is
  still accepted and means "same as label_min_fraction"). Argument checks are
  unclassed `stopifnot()`. The plot's default data is the per-node frame
  with `label_text` and `arc_fraction`, for introspection.
- **interactive**: `visHierarchicalLayout(direction = "LR", sortMethod =
  "directed", levelSeparation = 200, nodeSpacing = 120)`; dot nodes; cubic
  Bezier edges with arrows; `highlightNearest` on hover;
  `visInteraction(dragView = FALSE, zoomView = FALSE, dragNodes = TRUE)` —
  the widget is NOT zoomable (the old help page said it was). Node and edge
  sizes linearly rescaled into the pixel ranges (`.pt_scale_range()`); root
  colour `#222222`.
- The visNetwork check happens twice (dispatcher and backend), both
  `transitiontrees_missing_package`.
- `plot.transitiontrees_group()` adds `labs(caption = "group: <name>")` to
  ggplot results and `print()`s each plot as a side effect.
