# Fit a Prediction Suffix Tree from Categorical Sequence Data

Estimates a variable-depth context tree (prediction suffix tree; Ron,
Singer & Tishby 1996) from a collection of sequences. Each internal node
represents a context (string of recent states); each leaf carries a
smoothed conditional distribution over the next state. The tree is grown
to `max_depth`, then optionally pruned via
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md).

## Usage

``` r
context_tree(
  data,
  max_depth = 5L,
  min_count = 5L,
  smoothing = "floor",
  alphabet = NULL,
  weights = NULL,
  group = NULL,
  block = NULL,
  actor = NULL,
  time = NULL,
  action = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900
)
```

## Arguments

- data:

  Sequence data in any of these forms: a wide data.frame / character
  matrix (rows = trajectories, columns = time-steps), a list of
  character vectors, an `stslist`, or a **supported network/transition
  model object**, taken directly: a fitted network object carrying its
  sequences, or a transition-network object. Any such object that
  follows the same convention (a `$data`/`$sequences`/ `$seqdata`
  sequence slot) is detected structurally and works with no
  special-casing. For these objects the sequence frame is extracted from
  wherever the upstream stored it; an integer-coded frame is decoded
  through the object's label map (`$nodes` id/label table, positional
  `$labels`, or a `labels`/`alphabet` attribute), and that label set
  becomes the default alphabet so the tree shares the model's symbol
  space. A *pure graph* object that carries no sequences anywhere (an
  aggregated transition network) is rejected with guidance — sequences
  cannot be recovered from edge weights, the same reason numeric
  transition matrices are rejected; route to the sequence-bearing object
  explicitly.

- max_depth:

  Integer. Maximum context length the tree may represent. Default 5.

- min_count:

  Integer. Minimum number of times a context must occur to receive its
  own node. Default 5. Contexts seen fewer than `min_count` times are
  absorbed into their parent.

- smoothing:

  Smoothing specification: a method name as a string (uses defaults for
  that method's hyperparameters) or a list of the form
  `list(method, ...kwargs)` for explicit hyperparameters. Methods:
  `"floor"` (default; `ymin = 0.001`), `"laplace"` (`alpha = 1`),
  `"kneser_ney"` (`discount = 0.75`), `"witten_bell"`,
  `"jelinek_mercer"` (`lambda = 0.5`). The `"floor"` method also takes
  `rule`: `"interpolate"` (default, the interpolating floor — a
  distribution with a zero-count state is shifted toward uniform so each
  zero lands at exactly `ymin`) or `"cap"` (clamp every probability up
  to `ymin` and renormalise), e.g.
  `list("floor", ymin = 0.001, rule = "cap")`.

- alphabet:

  Character vector. Optional. Override the data-derived alphabet (useful
  when the test set may include states unseen in training).

- weights:

  Numeric vector of per-sequence weights, length equal to the number of
  input rows / list elements. If `NULL` (default) and `data` is an
  `stslist` carrying weights, those are auto-detected.

- group:

  Optional grouping for a **batch fit**: a vector with one entry per
  input sequence, a column name of a network object's `$metadata`, or —
  in long-format mode (`action` given) — a column name of `data`
  (collapsed to one value per sequence). When supplied (or when `data`
  is itself a grouped object such as a `netobject_group` or a
  `group_tna`), `context_tree()` fits one tree per group over a shared
  alphabet and returns a `transitiontrees_group` (a named list of
  `transitiontrees`s). Default `NULL` (single tree).

- block:

  Optional block id for a stratified group comparison, in the same
  shapes as `group` (per-sequence vector, or a column name of `data` in
  long-format mode). Carried on the returned `transitiontrees_group` and
  used automatically by
  [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
  so the caller never re-aligns it. Only meaningful alongside `group`.
  Default `NULL`.

- actor, time, order, session, time_threshold:

  Passed to
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  when `action` is given: the unit each sequence belongs to (`actor`),
  the timestamp column (`time`), an explicit ordering (`order`) or
  session id (`session`), and the gap in seconds that starts a new
  session (`time_threshold`, default 900).

- action:

  Character. Naming a state/code column switches `data` to
  **long-format** mode: it is reshaped to a wide sequence frame with
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  before fitting. Required to use any of `actor`/`time`/`order`/
  `session`. Default `NULL` (data is already in sequence shape).

## Value

For a single fit, a `transitiontrees` object (described below). For a
grouped fit (`group =` supplied, or a grouped family object passed in) a
`transitiontrees_group`: a named list of `transitiontrees`s, one per
group, in the group's key order, with its own `print` and
`as.data.frame` methods. A single `transitiontrees` is a list with
components

- nodes:

  Named list of node descriptors. Names are context strings (e.g.
  `"A -> B"`); the root is keyed by the literal sentinel `"<root>"`.
  Each entry has `depth` (integer), `counts` (numeric vector indexed by
  alphabet), `prob` (smoothed probability vector), and `n` (sum of
  counts).

- edges:

  data.frame with columns `parent`, `child`, `symbol` for fast tree
  traversal.

- alphabet:

  character vector of states.

- max_depth:

  Integer. The fitted depth (may be less than the requested `max_depth`
  if data are short).

- nmin:

  the chosen min-count threshold.

- smoothing:

  resolved smoothing list (`method` + hyperparameters).

- n_seq, n_obs:

  number of sequences and observations.

- pruned:

  Logical. `TRUE` after
  [`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
  has been applied.

- pruning:

  When `pruned` is `TRUE`, a list capturing the criterion / alpha /
  threshold used; otherwise `NULL`.

- data:

  The cleaned trajectories (a list of character vectors) retained for
  downstream bootstrap and permutation routines.

## Details

Construction is bottom-up via k-gram counting. The root holds the
marginal next-state distribution; depth-k nodes hold the next-state
distribution conditional on the most recent k states. Nodes whose total
count falls below `min_count` are not created. All nodes start "live"
and unpruned;
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
decides which to retain.

## References

Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia: learning
probabilistic automata with variable memory length. *Machine Learning*,
25, 117-149.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
tree
#> <transitiontrees>  40 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : A, B, C
#>   fit on   : 50 sequences, 600 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=600    -> A (0.37)
#> |-- A         n=200    -> A (0.38)
#> |   |-- A         n=66     -> A (0.38)
#> |   |   |-- A         n=18     -> A (0.50)
#> |   |   |-- B         n=31     -> C (0.39)
#> |   |   `-- C         n=10     -> A (0.50)
#> |   |-- B         n=67     -> A (0.48)
#> |   |   |-- A         n=22     -> A (0.64)
#> |   |   |-- B         n=17     -> A (0.41)
#> |   |   `-- C         n=21     -> A (0.38)
#> |   `-- C         n=49     -> C (0.41)
#> |       |-- A         n=15     -> C (0.40)
#> |       |-- B         n=15     -> C (0.53)
#> |       `-- C         n=17     -> B (0.41)
#> |-- B         n=181    -> A (0.40)
#> |   |-- A         n=56     -> A (0.43)
#> |   |   |-- A         n=20     -> B (0.60)
#> |   |   |-- B         n=15     -> A (0.33)
#> |   |   `-- C         n=16     -> A (0.56)
#> |   |-- B         n=58     -> B (0.43)
#> |   |   |-- A         n=20     -> B (0.70)
#> |   |   |-- B         n=22     -> A (0.41)
#> |   |   `-- C         n=12     -> B (0.42)
#> |   `-- C         n=51     -> A (0.43)
#> |       |-- A         n=18     -> A (0.39)
#> |       |-- B         n=8      -> A (0.62) 
#> ... 14 more nodes (use as.data.frame(x) or summary(x))
summary(tree)
#> <transitiontrees summary>  40 nodes, depth <= 3, 3 states  [unpruned]
#> 
#>  pathway depth count likely_next next_probability  divergence
#>  (start)     0   600           A        0.3666667          NA
#>        A     1   200           A        0.3800000 0.003825681
#>        B     1   181           A        0.3977901 0.013246404
#>        C     1   169           C        0.3491124 0.008133498
#>   B -> A     2    67           A        0.4776119 0.030761672
#>   A -> A     2    66           A        0.3787879 0.003552395
#>   A -> C     2    60           B        0.3500000 0.001283401
#>   B -> B     2    58           B        0.4310345 0.021097272
#>   A -> B     2    56           A        0.4285714 0.011114708
#>   C -> C     2    51           B        0.3725490 0.018405476
#>  changes_prediction
#>                  NA
#>               FALSE
#>               FALSE
#>                TRUE
#>               FALSE
#>               FALSE
#>                TRUE
#>                TRUE
#>               FALSE
#>                TRUE
#> # ... 30 more rows (use as.data.frame(tree) for the full table)
# }
```
