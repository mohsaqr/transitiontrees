# Changelog

## transitiontrees 0.1.2

CRAN release: 2026-06-18

Initial CRAN release.

### Fitting

- [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  fits a variable-depth pathway tree (prediction suffix tree; Ron,
  Singer & Tishby 1996) from a wide character matrix / data.frame, a
  list of character vectors, a long event log (`actor` / `time` /
  `action` / `order` / `session` arguments), an `stslist` state-sequence
  object, or a transition/network object.
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
  `"dendrogram"`, `"icicle"` (`ggraph`), and `"interactive"`
  (`visNetwork`).
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
