# Plot a Context Tree

Renders a fitted transitiontrees in one of four styles:

- `"horizontal"` (default) — pure-ggplot2 left-to-right phylogram: root
  on the left, contexts fanned out vertically to the right, each
  labelled beneath with its full arrow-notation context and (by default)
  its modal prediction `"(state pct%)"` on a second line. Node fill is
  the *most recent* move (rightmost token), so each branch off a depth-1
  hub shares a colour. Set `show_prediction = FALSE` for context-only
  labels.

- `"dendrogram"` — pure-ggplot2 radial tree: root at the centre, leaves
  on the outer ring.

- `"icicle"` — circular partition / sunburst via `ggraph` (Suggests);
  inner ring carries full state names, outer rings carry 3-letter
  abbreviations.

- `"interactive"` — `visNetwork` htmlwidget (Suggests). A draggable,
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

`"icicle"` requires `ggraph` + `tidygraph`; `"interactive"` requires
`visNetwork`. The dispatcher errors informatively if a needed Suggests
dependency is missing.

## Examples

``` r
# \donttest{
set.seed(1)
m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
tr <- context_tree(m, max_depth = 2L, min_count = 3L)
plot(tr)                           # left-to-right phylogram (default)

plot(tr, style = "dendrogram")     # radial dendrogram

if (requireNamespace("ggraph", quietly = TRUE) &&
    requireNamespace("tidygraph", quietly = TRUE))
  plot(tr, style = "icicle")       # sunburst (needs ggraph + tidygraph)

if (requireNamespace("visNetwork", quietly = TRUE))
  plot(tr, style = "interactive")  # visNetwork htmlwidget (Suggests)

{"x":{"nodes":{"id":["<root>","A","C","B","A -> C","C -> B","B -> A","A -> A","A -> B","B -> B","C -> A","C -> C","B -> C"],"label":["(start)","A","C","B","A","C","B","A","A","B","C","C","B"],"count":[200,58,59,63,22,22,24,13,17,17,14,15,16],"level":[0,1,1,1,2,2,2,2,2,2,2,2,2],"color":["#222222","#ED90A4","#6FB1E7","#7EBA68","#ED90A4","#6FB1E7","#7EBA68","#ED90A4","#ED90A4","#7EBA68","#6FB1E7","#6FB1E7","#7EBA68"],"title":["<b>Pathway<\/b>: (start)<br/><b>Depth<\/b>: 0<br/><b>Count<\/b>: 200<br/><b>Modal next<\/b>: B (0.340)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.330<br/>&nbsp;&nbsp;B: 0.340<br/>&nbsp;&nbsp;C: 0.330","<b>Pathway<\/b>: A<br/><b>Depth<\/b>: 1<br/><b>Count<\/b>: 58<br/><b>Modal next<\/b>: C (0.431)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.259<br/>&nbsp;&nbsp;B: 0.310<br/>&nbsp;&nbsp;C: 0.431","<b>Pathway<\/b>: C<br/><b>Depth<\/b>: 1<br/><b>Count<\/b>: 59<br/><b>Modal next<\/b>: B (0.424)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.288<br/>&nbsp;&nbsp;B: 0.424<br/>&nbsp;&nbsp;C: 0.288","<b>Pathway<\/b>: B<br/><b>Depth<\/b>: 1<br/><b>Count<\/b>: 63<br/><b>Modal next<\/b>: A (0.429)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.429<br/>&nbsp;&nbsp;B: 0.286<br/>&nbsp;&nbsp;C: 0.286","<b>Pathway<\/b>: A -> C<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 22<br/><b>Modal next<\/b>: B (0.409)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.273<br/>&nbsp;&nbsp;B: 0.409<br/>&nbsp;&nbsp;C: 0.318","<b>Pathway<\/b>: C -> B<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 22<br/><b>Modal next<\/b>: A (0.500)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.500<br/>&nbsp;&nbsp;B: 0.136<br/>&nbsp;&nbsp;C: 0.364","<b>Pathway<\/b>: B -> A<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 24<br/><b>Modal next<\/b>: C (0.458)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.250<br/>&nbsp;&nbsp;B: 0.292<br/>&nbsp;&nbsp;C: 0.458","<b>Pathway<\/b>: A -> A<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 13<br/><b>Modal next<\/b>: C (0.538)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.154<br/>&nbsp;&nbsp;B: 0.308<br/>&nbsp;&nbsp;C: 0.538","<b>Pathway<\/b>: A -> B<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 17<br/><b>Modal next<\/b>: A (0.529)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.529<br/>&nbsp;&nbsp;B: 0.235<br/>&nbsp;&nbsp;C: 0.235","<b>Pathway<\/b>: B -> B<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 17<br/><b>Modal next<\/b>: B (0.412)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.294<br/>&nbsp;&nbsp;B: 0.412<br/>&nbsp;&nbsp;C: 0.294","<b>Pathway<\/b>: C -> A<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 14<br/><b>Modal next<\/b>: B (0.500)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.214<br/>&nbsp;&nbsp;B: 0.500<br/>&nbsp;&nbsp;C: 0.286","<b>Pathway<\/b>: C -> C<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 15<br/><b>Modal next<\/b>: B (0.467)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.400<br/>&nbsp;&nbsp;B: 0.467<br/>&nbsp;&nbsp;C: 0.133","<b>Pathway<\/b>: B -> C<br/><b>Depth<\/b>: 2<br/><b>Count<\/b>: 16<br/><b>Modal next<\/b>: C (0.500)<br/><b>Next-state distribution<\/b>:<br/>&nbsp;&nbsp;A: 0.188<br/>&nbsp;&nbsp;B: 0.312<br/>&nbsp;&nbsp;C: 0.500"],"size":[45,18.42245989304813,18.6096256684492,19.35828877005348,11.68449197860962,11.68449197860962,12.05882352941176,10,10.74866310160428,10.74866310160428,10.18716577540107,10.37433155080214,10.56149732620321]},"edges":{"from":["<root>","<root>","<root>","C","B","A","A","B","B","A","C","C"],"to":["A","C","B","A -> C","C -> B","B -> A","A -> A","A -> B","B -> B","C -> A","C -> C","B -> C"],"width":[9.1,9.280000000000001,10,2.62,2.62,2.98,1,1.72,1.72,1.18,1.36,1.54]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","borderWidth":1,"color":{"border":"#444444","highlight":{"border":"#000000"}},"font":{"size":16,"face":"sans-serif"}},"manipulation":{"enabled":false},"layout":{"hierarchical":{"enabled":true,"levelSeparation":200,"nodeSpacing":120,"direction":"LR","sortMethod":"directed"}},"edges":{"arrows":"to","color":{"color":"#cfcfcf","highlight":"#555555"},"smooth":{"enabled":true,"type":"cubicBezier","roundness":0.5}},"interaction":{"dragNodes":true,"dragView":false,"hover":true,"zoomView":false,"zoomSpeed":1}},"groups":null,"width":null,"height":null,"idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false},"main":null,"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":true,"hoverNearest":true,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}# }
```
