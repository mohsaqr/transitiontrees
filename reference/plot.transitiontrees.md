# Plot a Context Tree

Renders a fitted transitiontrees in one of four styles:

- `"horizontal"` (default) – pure-ggplot2 left-to-right phylogram: root
  on the left, contexts fanned out vertically to the right, each
  labelled beneath with its full arrow-notation context and (by default)
  its modal prediction `"(state pct%)"` on a second line. Node fill is
  the *most recent* move (rightmost token), so each branch off a depth-1
  hub shares a colour. Set `show_prediction = FALSE` for context-only
  labels.

- `"dendrogram"` – pure-ggplot2 radial tree: root at the centre, leaves
  on the outer ring.

- `"icicle"` – pure-ggplot2 circular partition / sunburst; inner ring
  carries full state names, outer rings carry 3-letter abbreviations.

- `"interactive"` – `visNetwork` htmlwidget (Suggests). A draggable,
  zoomable hierarchical tree; node size = context count and edge width =
  child's count (“flow”), the same encoding as the static styles. Hover
  for a tooltip with the full pathway, count, modal next state, and the
  complete next-state distribution.

Common encoding: node size = context count, edge thickness = child's
count (“flow”). Node fill is the most recent move (rightmost token of
the context) in the `"horizontal"` style; the `"dendrogram"`,
`"icicle"`, and `"interactive"` styles colour by the branching (oldest)
token.

## Usage

``` r
# S3 method for class 'transitiontrees'
plot(
  x,
  style = c("horizontal", "dendrogram", "icicle", "interactive"),
  point_size_range = NULL,
  edge_size_range = NULL,
  ...
)
```

## Arguments

- x:

  A `transitiontrees`.

- style:

  One of `"horizontal"` (default), `"dendrogram"`, `"icicle"`, or
  `"interactive"`.

- point_size_range:

  Numeric length-2 vector controlling the minimum and maximum node-point
  size. Default `c(5, 16)` for `"dendrogram"`, `c(4, 14)` for
  `"horizontal"`, `c(10, 45)` (pixels) for `"interactive"`. Ignored by
  `"icicle"`.

- edge_size_range:

  Numeric length-2 vector for edge width. Default `c(0.3, 2.5)` for the
  static styles, `c(1, 10)` (pixels) for `"interactive"`. Ignored by
  `"icicle"`.

- ...:

  Passed to the chosen backend. For `style = "horizontal"`,
  `show_prediction` (logical, default `TRUE`) toggles each node's modal
  prediction `"(state pct%)"` on a second label line; set `FALSE` for
  context-only labels. For `style = "interactive"`, `width` / `height`
  size the htmlwidget.

## Value

A ggplot object for the three static styles; an `htmlwidget` for
`"interactive"`.

## Details

The three static styles (`"horizontal"`, `"dendrogram"`, `"icicle"`) are
drawn in pure ggplot2. `"interactive"` requires `visNetwork`; the
dispatcher errors informatively if it is missing.

## Examples

``` r
set.seed(1)
m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr <- context_tree(m, max_depth = 2L, min_count = 3L)
plot(tr)                           # left-to-right phylogram (default)

plot(tr, style = "dendrogram")     # radial dendrogram

plot(tr, style = "icicle")         # sunburst

if (requireNamespace("visNetwork", quietly = TRUE)) {
  widget <- plot(tr, style = "interactive")  # visNetwork htmlwidget
  if (interactive()) print(widget)           # opens in the viewer
}
```
