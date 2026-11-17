# Implementation notes removed from `?plot_trajectories` (2026-11-17)

Still true of R/plot_trajectories.R; kept for maintainers.

- `.trajectory_data()` counts every prefix of every sequence (`table()`),
  keeps those with count >= `min_count`, drops any prefix whose parent was
  dropped (so the tree stays connected). The root count is the number of
  sequences.
- `eprob`: for a first move, the root distribution; otherwise
  `query_pathway(tree, <last max_depth states of the parent prefix>)`, which
  falls back to the longest fitted suffix. A state outside the alphabet
  yields `NA`.
- Layout: leaves stacked at y = 1, 2, ...; parents at the mean of their
  children; 40-point cosine-smoothed edges.
- Node glyphs: capsules traced by hand (half-extents 0.17 x 0.26, corner
  radius 0.17, 6 points per corner); text white when the rescaled value
  > 0.45.
- Colour: `scale_fill_viridis_c(option = "C", direction = -1)` (higher =
  darker); predictability limits fixed to `[0, 1]`. Edge width range
  `c(0.4, 3.6)`.
