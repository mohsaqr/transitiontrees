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
  character vectors, or a **supported network/transition model object**,
  taken directly: a fitted network object carrying its sequences, or a
  transition-network object. Any such object that follows the same
  convention (a `$data`/`$sequences`/ `$seqdata` sequence slot) is
  detected structurally and works with no special-casing. For these
  objects the sequence frame is extracted from wherever the upstream
  stored it; an integer-coded frame is decoded through the object's
  label map (`$nodes` id/label table, positional `$labels`, or a
  `labels`/`alphabet` attribute), and that label set becomes the default
  alphabet so the tree shares the model's symbol space. A *pure graph*
  object that carries no sequences anywhere (an aggregated transition
  network) is rejected with guidance – sequences cannot be recovered
  from edge weights, the same reason numeric transition matrices are
  rejected; route to the sequence-bearing object explicitly. In wide /
  list input, `NA`, the empty string, and the TraMineR codes `"%"`
  (void) and `"*"` (missing) all mark "no observation": the cell is
  dropped (gaps closed, trailing padding trimmed) rather than treated as
  a state.

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
  `rule`: `"interpolate"` (default, the interpolating floor – a
  distribution with a zero-count state is shifted toward uniform so each
  zero lands at exactly `ymin`) or `"cap"` (clamp every probability up
  to `ymin` and renormalise), e.g.
  `list("floor", ymin = 0.001, rule = "cap")`.

- alphabet:

  Character vector. Optional. Override the data-derived alphabet (useful
  when the test set may include states unseen in training).

- weights:

  Numeric vector of per-sequence weights, length equal to the number of
  input rows / list elements. `NULL` (default) fits unweighted.

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

  the `min_count` threshold used for the fit.

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

- weights:

  The per-sequence weights aligned to `data`, or `NULL` for an
  unweighted fit. Weighted trees are refused by
  [`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
  [`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
  and
  [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
  which rebuild counts from the raw sequences.

Read the fitted tree through its accessors rather than these components:
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
/ [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for the
node table,
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
for the fit scalars, and
[`n_nodes()`](https://pak.dynasite.org/transitiontrees/reference/n_nodes.md)
for the size.

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
[doi:10.1023/A:1026490906255](https://doi.org/10.1023/A%3A1026490906255)

## Examples

``` r
# Wide character matrix: one learner per row, one time-step per column
tree <- context_tree(trajectories, max_depth = 3L, min_count = 5L)
tree
#> <transitiontrees>  38 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |-- Average   n=70     -> Active (0.53)
#> |   |   `-- Disengaged  n=15     -> Active (0.67)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   |   |-- Active    n=53     -> Average (0.55)
#> |   |   |-- Average   n=66     -> Active (0.58)
#> |   |   `-- Disengaged  n=12     -> Average (0.83)
#> |   `-- Disengaged  n=29     -> Active (0.52)
#> |       |-- Active    n=6      -> Active (0.83)
#> |       |-- Average   n=13     -> Active (0.54)
#> |       `-- Disengaged  n=9      -> Average (0.56)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   |-- Active    n=74     -> Average (0.55)
#> |   |   |-- Average   n=64     -> Active (0.45)
#> |   |   `-- Disengaged  n=10     -> Average (0.50)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.57)
#> |   |   |-- Average   n=248    -> Average (0.72)
#> |   |   `-- Disengaged  n=63     -> Average (0.65)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> |       |-- Active    n=7      -> Disengaged (0.57)
#> |       |-- Average   n=64     -> Disengaged (0.48) 
#> ... 12 more nodes (use as.data.frame(x) or summary(x))

# A list of character vectors works the same way
set.seed(1)
seqs <- replicate(50, sample(c("A", "B", "C"), 12, replace = TRUE),
                  simplify = FALSE)
context_tree(seqs, max_depth = 2L)
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : A, B, C
#>   fit on   : 50 sequences, 600 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=600    -> A (0.37)
#> |-- A         n=200    -> A (0.38)
#> |   |-- A         n=66     -> A (0.38)
#> |   |-- B         n=67     -> A (0.48)
#> |   `-- C         n=49     -> C (0.41)
#> |-- B         n=181    -> A (0.40)
#> |   |-- A         n=56     -> A (0.43)
#> |   |-- B         n=58     -> B (0.43)
#> |   `-- C         n=51     -> A (0.43)
#> `-- C         n=169    -> C (0.35)
#>     |-- A         n=60     -> B (0.35)
#>     |-- B         n=42     -> C (0.43)
#>     `-- C         n=51     -> B (0.37) 

# Long event log: name the columns and it is reshaped first
context_tree(group_regulation_long, actor = "Actor", time = "Time",
             action = "Action", max_depth = 2L)
#> <transitiontrees>  87 nodes, depth <= 2, 9 states  [unpruned]
#>   alphabet : adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis
#>   fit on   : 2000 sequences, 27533 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=27533  -> consensus (0.25)
#> |-- adapt     n=509    -> consensus (0.47)
#> |   |-- consensus  n=27     -> cohesion (0.48)
#> |   |-- coregulate  n=28     -> consensus (0.50)
#> |   |-- discuss   n=259    -> consensus (0.47)
#> |   |-- emotion   n=6      -> consensus (0.83)
#> |   |-- monitor   n=16     -> consensus (0.50)
#> |   |-- plan      n=6      -> consensus (0.66)
#> |   `-- synthesis  n=140    -> consensus (0.48)
#> |-- cohesion  n=1695   -> consensus (0.50)
#> |   |-- adapt     n=130    -> consensus (0.54)
#> |   |-- cohesion  n=45     -> consensus (0.42)
#> |   |-- consensus  n=84     -> consensus (0.51)
#> |   |-- coregulate  n=61     -> consensus (0.50)
#> |   |-- discuss   n=173    -> consensus (0.54)
#> |   |-- emotion   n=843    -> consensus (0.47)
#> |   |-- monitor   n=71     -> consensus (0.53)
#> |   |-- plan      n=147    -> consensus (0.52)
#> |   `-- synthesis  n=20     -> consensus (0.55)
#> |-- consensus  n=6329   -> plan (0.40)
#> |   |-- adapt     n=222    -> plan (0.37)
#> |   |-- cohesion  n=795    -> plan (0.38)
#> |   |-- consensus  n=483    -> plan (0.39)
#> |   |-- coregulate  n=252    -> plan (0.40)
#> |   |-- discuss   n=1180   -> plan (0.41)
#> |   |-- emotion   n=830    -> plan (0.39) 
#> ... 61 more nodes (use as.data.frame(x) or summary(x))

# One tree per group over a shared alphabet
context_tree(trajectories, group = rep(c("early", "late"), each = 69),
             max_depth = 2L)
#> <transitiontrees_group>  2 groups
#>   early           13 nodes, depth <= 2, 69 seq, 969 obs
#>   late            13 nodes, depth <= 2, 67 seq, 901 obs
```
