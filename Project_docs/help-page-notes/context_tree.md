# Implementation notes removed from `?context_tree` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/context_tree.R, R/prepare_input.R); it is kept
here for maintainers.

- **Network objects** (`.ct_is_netobject()`, `.ct_unwrap_netobject()`):
  detected by known class, or structurally: any list-like object with a
  2-D sequence frame in a `$data`, `$sequences` or `$seqdata` slot (also an
  embedded netobject, or a `"data"`/`"sequences"` attribute). The `$data`
  slot is tried first.
- **Integer-coded sequence frames** inside a network object are decoded,
  not rejected: codes are mapped through a `$nodes` id/label table
  (`label[match(code, id)]`), else a positional `$labels` vector, else a
  `labels`/`alphabet` attribute, else cast to character. `NA` stays
  end-of-sequence. The decoded label set becomes the default alphabet.
  The "reject numeric matrices" rule applies only to a top-level numeric
  matrix (which collides with "square numeric = transition matrix" in the
  sibling network packages).
- **Pure graph objects** (nodes/edges/weights with the sequences absent
  everywhere, e.g. an aggregated transition network whose constructor
  nulls the data slot) raise `transitiontrees_bad_input` with guidance:
  sequences cannot be recovered from edge weights.
- **Missing codes** (`.CT_MISSING`): `NA`, `""`, TraMineR `"%"` (void)
  and `"*"` (missing). Internal gaps are closed and trailing padding is
  trimmed; the same cleaning rule is used for list and matrix input, so
  identical data gives identical trees whatever the input shape.
- **Grouped fits** (`.ct_fit_grouped()`): a grouped family object
  (`netobject_group`, `group_tna`, or a named list whose every element is
  a fittable network object) dispatches over its elements; `group =`
  splits one dataset. Passing both errors. `weights` with a grouped
  object errors. A single-string `group` names a `$metadata` column for a
  network object, or a column of `data` in long mode. The alphabet is the
  explicit `alphabet` or the union of the groups' symbol sets. The
  `transitiontrees_group` has class `c("transitiontrees_group", "list")`
  and a `group` attribute (the grouping-variable name or `NA`); members
  carry an aligned `block` attribute when `block` is given.
- **Long mode**: `group`/`block` given as column names are carried
  through `prepare_input(meta = ...)` (first event's value per session),
  so they are aligned to the wide rows by construction, for every
  reshape path (time, session, order).
- **Weights**: validated and reindexed to the trajectories that survive
  cleaning; uniform weights are treated as unweighted (`NULL`); at least
  one positive weight is required. Weighted trees are refused by
  `bootstrap_pathways()`, `compare_trees()` and `compare_groups()`
  (`transitiontrees_weighted_unsupported`), which rebuild counts from the
  raw sequences.
- **Construction**: k-gram counting top-down; each depth is smoothed
  after its parent so back-off smoothers can read the parent's
  probabilities. If no context of a depth clears `min_count`, the fitted
  `max_depth` is lowered.
- **Object components** (`transitiontrees`): `nodes` (named list keyed by
  context string; the root key is the internal sentinel `"<root>"`,
  displayed as `"(start)"`; each node has `depth`, `counts`, `prob`, `n`),
  `edges` (data.frame `parent`, `child`, `symbol`), `alphabet`,
  `max_depth` (fitted), `nmin` (the `min_count` used), `smoothing`
  (resolved list), `n_seq`, `n_obs` (weighted token total = sum of root
  counts), `pruned`, `pruning` (criterion/alpha/threshold or `NULL`),
  `data` (cleaned trajectories, kept for bootstrap/permutation), `weights`.
- **PST parity**: with `smoothing = list("floor", ymin = ..., rule =
  "interpolate")` (the default), counts and node probabilities match PST
  0.94.1 at all depths (see PARITY.md). `rule = "cap"` is the package's
  original clamp-and-renormalise rule.
