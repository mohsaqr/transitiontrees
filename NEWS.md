# transitiontrees 0.2.0

This release collects the changes made since the CRAN release 0.1.2
(including the unreleased 0.1.3 and 0.1.4).

## New functions and arguments

* `state_frequencies()` returns the count and proportion of each state
  (one table per group for a grouped fit).
* `summary()` methods for groups of trees, two-tree comparisons and tuning
  grids; `summary()` of a group comparison now returns a one-row summary
  object with its own print method instead of printing and returning the
  pathway table invisibly.
* `as.data.frame()` methods gain a `what` argument to reach every table of a
  result: `what = "null"` for the permutation null of `compare_trees()`,
  `what = "omnibus"` and `what = "distance"` for `compare_groups()`,
  `what = "best"` for `tune_tree()` (new method), and `what = "resamples"`
  (with `stat`) for `bootstrap_pathways()`.
* `predict()` gains `tidy = TRUE`, returning one row per history with the
  matched context, the most likely next state and the full distribution.
  `generate_sequences()` and `simulate()` gain `tidy = TRUE`, returning a
  wide data frame that `context_tree()` accepts. The default outputs are
  unchanged.
* `context_tree()` accepts the object returned by `tna::prepare_data()`.

## Changes in behavior

* The column of the `tune_tree()` grid that holds the minimum count is now
  named `min_count`, matching the argument (it was `nmin`).
* A tuning configuration whose folds all fail now reports `logLik`,
  `perplexity` and `n_nodes_avg` as `NA` (the log-likelihood was 0). The
  warning includes the first fold error.
* `predict()` splits a history given as an arrow-notation string
  (`"A -> B"`) into its states, and warns when a history contains a state
  outside the alphabet.
* Errors and warnings carry condition classes (`transitiontrees_bad_input`,
  `transitiontrees_alphabet_mismatch`, `transitiontrees_no_data`,
  `transitiontrees_weighted_unsupported`, `transitiontrees_empty_result`,
  `transitiontrees_missing_package`, and the warnings
  `transitiontrees_fold_failed`, `transitiontrees_no_best`,
  `transitiontrees_unknown_state`,
  `transitiontrees_degenerate_permutation`).
* `plot(tune_tree(...))` marks the selected configuration only in its own
  facet; `plot()` of a bootstrap result no longer emits a ggplot2 message.
* `plot_pruning()` labels use ASCII text, so the plot draws in any locale.
* ggplot2 (>= 3.5.0) is required.

## Documentation

* New vignettes "Variable-order transition trees" and "Comparing groups".
  Eight guides and three case studies on the package website replace the
  previous five vignettes.

# transitiontrees 0.1.2

Initial CRAN release.

## Fitting

* `context_tree()` fits a variable-depth pathway tree (prediction
  suffix tree; Ron, Singer & Tishby 1996) from a wide character
  matrix / data.frame, a list of character vectors, a long event log
  (`actor` / `time` / `action` / `order` / `session` arguments), or a
  transition/network object. `NA`, `""`, and the TraMineR codes `%`
  (void) and `*` (missing) are treated as gaps in wide / list input.
* `prepare_input()` reshapes a long event log to a wide sequence frame
  (timestamp / session logic), and can carry per-sequence metadata
  through the reshape via `meta`.
* Five smoothing schemes via the unified `smoothing` argument
  (`"floor"`, `"laplace"`, `"kneser_ney"`, `"witten_bell"`,
  `"jelinek_mercer"`); hyperparameters as `list(method, ...)`.
* `prune_tree()` supports four criteria: likelihood-ratio `G2`,
  Kullback-Leibler, `AIC`, `BIC`.
* `smooth_tree()` re-smooths a fitted tree; `model_fit()` /
  `n_nodes()` are tidy fit-summary accessors.
* Grouped fits: `context_tree(..., group =)` (a per-sequence vector or
  a column name) fits one tree per group over a shared alphabet and
  returns a `transitiontrees_group`. `block =` carries a stratifying id
  (e.g. subject) for `compare_groups()`.

## Pathway-centric API

* `tree_pathways()`, `common_pathways()`, `divergent_pathways()`,
  `sharp_pathways()` rank pathways by frequency, divergence from the
  suffix-parent, or modal-flip status.
* `tree_dependence()` is the per-context entropy/divergence diagnostic
  table; `query_pathway()`, `pathway_exists()`, `subtree()` provide
  tree introspection.

## Prediction, scoring, and imputation

* `predict()` / `simulate()` / `generate_sequences()` for next-state
  prediction and sampling.
* `logLik()`, `nobs()`, `AIC()`, `BIC()`, `perplexity()`,
  `score_sequences()`, `score_positions()` form the predictive-
  evaluation toolchain.
* `impute_sequences()` fills internal gaps in incomplete sequences.
* `mine_contexts()` / `mine_sequences()` scan for contexts where a
  state is unusually likely or unlikely and for the best/worst-fit held-out
  sequences.

## Resampling and group comparison

* `tune_tree()` k-fold cross-validates `max_depth`, `min_count`,
  smoothing, and pruning.
* `bootstrap_pathways()` reports per-pathway stability and
  informativeness with bootstrap CIs.
* `compare_trees()` runs a permutation test for two-tree divergence.
* `compare_groups()` compares a `transitiontrees_group` on two axes ---
  behavioral (Jensen-Shannon divergence of next-state distributions)
  and usage (prevalence) --- with a permutation null (optionally
  stratified by `block` for repeated-measures designs), Benjamini-
  Hochberg FDR, and a between-group distance matrix.
* `tree_distance()` computes count-weighted symmetric KL between two
  trees.

## Visualisation

* `plot()` on a `transitiontrees` offers four styles: `"horizontal"`
  (default), `"dendrogram"`, and `"icicle"` (all pure ggplot2), plus
  `"interactive"` (`visNetwork`). `plot()` on a `transitiontrees_group`
  draws one figure per group.
* `plot_pathways()`, `plot_divergence()`, `plot_distributions()`,
  `plot_predictive()`, `plot_pathway_resamples()`, and the
  bootstrap / comparison / tuning plot methods.
* `plot_difference()` renders the early-vs-late style difference
  between two groups as a per-context map (Pearson residuals against
  the no-difference null, or raw probability difference) or on the
  context-tree layout.

## Bundled data

* `trajectories`, `group_regulation_long`, `ai_long`, and
  `engagement` for examples and tests.

## Validation

* Equivalence-tested at machine precision against an independent
  external reference implementation of the model --- counts exact,
  probabilities within 1.11e-16 --- and cross-checked against an
  independent first-order Markov reference and a standard long-to-sequence
  reshaper. The equivalence suite lives outside the package and is run
  locally.
