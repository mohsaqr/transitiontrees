# transitiontrees 0.1.3

## Documentation

* Every help page now has examples that run unwrapped (no `\donttest{}`),
  with fixed seeds, and examples were added for the `print()`, `summary()`,
  `plot()` and `as.data.frame()` methods, `simulate()` and `nobs()`.
* `@return` sections describe the returned structure (one row per what,
  every column), including the bootstrap summary columns, the
  `compare_groups()` omnibus table and `stratified` flag, and the
  `weights` component of a fitted tree.
* Corrected stale or wrong statements: `compare_pruning(alpha =)` is used
  by `"G2"` only (not by AIC/BIC); `compare_smoothing()` holds
  `min_count` (not `nmin`) fixed; the bootstrap G-squared level is the
  `alpha` argument; the `trajectories` help no longer points at a
  non-existent vignette; `plot_pruning()` no longer cites the unshipped
  `PARITY.md`.
* Added method references (AIC, BIC, smoothing, Jensen-Shannon divergence,
  Benjamini-Hochberg, bootstrap) and replaced non-ASCII punctuation in the
  help pages.
* README and vignettes: examples run on the bundled data, and results are
  read through the package's printed objects and accessors instead of
  `$` / `attr()` extraction.

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
