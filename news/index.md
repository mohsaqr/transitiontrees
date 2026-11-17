# Changelog

## transitiontrees 0.2.0

This release collects the changes made since the CRAN release 0.1.2
(including the unreleased 0.1.3 and 0.1.4).

### New functions and arguments

- [`state_frequencies()`](https://pak.dynasite.org/transitiontrees/reference/state_frequencies.md)
  returns the count and proportion of each state (one table per group
  for a grouped fit).
- [`summary()`](https://rdrr.io/r/base/summary.html) methods for groups
  of trees, two-tree comparisons and tuning grids;
  [`summary()`](https://rdrr.io/r/base/summary.html) of a group
  comparison now returns a one-row summary object with its own print
  method instead of printing and returning the pathway table invisibly.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  gain a `what` argument to reach every table of a result:
  `what = "null"` for the permutation null of
  [`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
  `what = "omnibus"` and `what = "distance"` for
  [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
  `what = "best"` for
  [`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
  (new method), and `what = "resamples"` (with `stat`) for
  [`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).
- [`predict()`](https://rdrr.io/r/stats/predict.html) gains
  `tidy = TRUE`, returning one row per history with the matched context,
  the most likely next state and the full distribution.
  [`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
  and [`simulate()`](https://rdrr.io/r/stats/simulate.html) gain
  `tidy = TRUE`, returning a wide data frame that
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  accepts. The default outputs are unchanged.
- [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  accepts the object returned by
  [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html).

### Changes in behavior

- The column of the
  [`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
  grid that holds the minimum count is now named `min_count`, matching
  the argument (it was `nmin`).
- A tuning configuration whose folds all fail now reports `logLik`,
  `perplexity` and `n_nodes_avg` as `NA` (the log-likelihood was 0). The
  warning includes the first fold error.
- [`predict()`](https://rdrr.io/r/stats/predict.html) splits a history
  given as an arrow-notation string (`"A -> B"`) into its states, and
  warns when a history contains a state outside the alphabet.
- Errors and warnings carry condition classes
  (`transitiontrees_bad_input`, `transitiontrees_alphabet_mismatch`,
  `transitiontrees_no_data`, `transitiontrees_weighted_unsupported`,
  `transitiontrees_empty_result`, `transitiontrees_missing_package`, and
  the warnings `transitiontrees_fold_failed`, `transitiontrees_no_best`,
  `transitiontrees_unknown_state`,
  `transitiontrees_degenerate_permutation`).
- `plot(tune_tree(...))` marks the selected configuration only in its
  own facet; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  a bootstrap result no longer emits a ggplot2 message.
- [`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md)
  labels use ASCII text, so the plot draws in any locale.
- ggplot2 (\>= 3.5.0) is required.

### Documentation

- New vignettes “Variable-order transition trees” and “Comparing
  groups”. Eight guides and three case studies on the package website
  replace the previous five vignettes.

## transitiontrees 0.1.2

CRAN release: 2026-06-18

Initial CRAN release.

### Fitting

- [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  fits a variable-depth pathway tree (prediction suffix tree; Ron,
  Singer & Tishby 1996) from a wide character matrix / data.frame, a
  list of character vectors, a long event log (`actor` / `time` /
  `action` / `order` / `session` arguments), or a transition/network
  object. `NA`, `""`, and the TraMineR codes `%` (void) and `*`
  (missing) are treated as gaps in wide / list input.
- [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  reshapes a long event log to a wide sequence frame (timestamp /
  session logic), and can carry per-sequence metadata through the
  reshape via `meta`.
- Five smoothing schemes via the unified `smoothing` argument
  (`"floor"`, `"laplace"`, `"kneser_ney"`, `"witten_bell"`,
  `"jelinek_mercer"`); hyperparameters as `list(method, ...)`.
- [`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
  supports four criteria: likelihood-ratio `G2`, Kullback-Leibler,
  `AIC`, `BIC`.
- [`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
  re-smooths a fitted tree;
  [`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
  /
  [`n_nodes()`](https://pak.dynasite.org/transitiontrees/reference/n_nodes.md)
  are tidy fit-summary accessors.
- Grouped fits: `context_tree(..., group =)` (a per-sequence vector or a
  column name) fits one tree per group over a shared alphabet and
  returns a `transitiontrees_group`. `block =` carries a stratifying id
  (e.g. subject) for
  [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

### Pathway-centric API

- [`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
  [`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
  [`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md),
  [`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md)
  rank pathways by frequency, divergence from the suffix-parent, or
  modal-flip status.
- [`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md)
  is the per-context entropy/divergence diagnostic table;
  [`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
  [`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md),
  [`subtree()`](https://pak.dynasite.org/transitiontrees/reference/subtree.md)
  provide tree introspection.

### Prediction, scoring, and imputation

- [`predict()`](https://rdrr.io/r/stats/predict.html) /
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) /
  [`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
  for next-state prediction and sampling.
- [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
  [`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
  [`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
  form the predictive- evaluation toolchain.
- [`impute_sequences()`](https://pak.dynasite.org/transitiontrees/reference/impute_sequences.md)
  fills internal gaps in incomplete sequences.
- [`mine_contexts()`](https://pak.dynasite.org/transitiontrees/reference/mine_contexts.md)
  /
  [`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md)
  scan for contexts where a state is unusually likely or unlikely and
  for the best/worst-fit held-out sequences.

### Resampling and group comparison

- [`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
  k-fold cross-validates `max_depth`, `min_count`, smoothing, and
  pruning.
- [`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
  reports per-pathway stability and informativeness with bootstrap CIs.
- [`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
  runs a permutation test for two-tree divergence.
- [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
  compares a `transitiontrees_group` on two axes — behavioral
  (Jensen-Shannon divergence of next-state distributions) and usage
  (prevalence) — with a permutation null (optionally stratified by
  `block` for repeated-measures designs), Benjamini- Hochberg FDR, and a
  between-group distance matrix.
- [`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)
  computes count-weighted symmetric KL between two trees.

### Visualisation

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `transitiontrees` offers four styles: `"horizontal"` (default),
  `"dendrogram"`, and `"icicle"` (all pure ggplot2), plus
  `"interactive"` (`visNetwork`).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `transitiontrees_group` draws one figure per group.
- [`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
  [`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md),
  [`plot_distributions()`](https://pak.dynasite.org/transitiontrees/reference/plot_distributions.md),
  [`plot_predictive()`](https://pak.dynasite.org/transitiontrees/reference/plot_predictive.md),
  [`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md),
  and the bootstrap / comparison / tuning plot methods.
- [`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
  renders the early-vs-late style difference between two groups as a
  per-context map (Pearson residuals against the no-difference null, or
  raw probability difference) or on the context-tree layout.

### Bundled data

- `trajectories`, `group_regulation_long`, `ai_long`, and `engagement`
  for examples and tests.

### Validation

- Equivalence-tested at machine precision against an independent
  external reference implementation of the model — counts exact,
  probabilities within 1.11e-16 — and cross-checked against an
  independent first-order Markov reference and a standard
  long-to-sequence reshaper. The equivalence suite lives outside the
  package and is run locally.
