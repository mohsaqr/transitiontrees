# transitiontrees

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/transitiontrees)](https://CRAN.R-project.org/package=transitiontrees)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/transitiontrees)](https://cran.r-project.org/package=transitiontrees)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Predictive pathway discovery in categorical sequence data.

`transitiontrees` fits a variable-depth pathway tree (a prediction
suffix tree; Ron, Singer & Tishby 1996) from sequences and exposes a
tidy, pathway-centric API: most common pathways, most predictively
divergent pathways, modal-flip diagnostics, and next-state predictions.

## Why this package

Sequence-analysis tooling tends to either count k-grams (descriptive)
or fit fixed-order Markov chains (one-size-fits-all memory). Real
trajectories — student engagement weeks, patient pathways, clickstreams,
play sequences — usually have memory that varies *by context*: some
states predict the next state on their own, others only when paired
with a longer history. `transitiontrees` makes that variable-depth
structure the central object and reports it as a ranked list of pathways
the data actually supports.

It is a self-contained, independent implementation of the variable-order
prediction-suffix-tree model:

- pure base R with a single plotting dependency (ggplot2); no compiled
  code and no heavy framework dependencies;
- tidy `data.frame` outputs by default, in one canonical schema so every
  table joins cleanly;
- validated at machine precision against independent external reference
  implementations of the model (counts, node probabilities at all depths,
  context queries, per-position prediction, `logLik`, and topology), as
  documented in
  [`PARITY.md`](https://github.com/mohsaqr/transitiontrees/blob/main/PARITY.md)
  (the equivalence suite runs locally and is not shipped with the package);
- a pathway-centric API (`tree_pathways()`, `common_pathways()`,
  `divergent_pathways()`, `sharp_pathways()`) that ranks trajectories
  by frequency, predictive divergence, or modal-flip — the structure
  domain experts actually want to read off the model;
- verbs deliberately namespaced (`prune_tree()`, `tree_pathways()`,
  `tree_dependence()`) so `library(transitiontrees)` does not shadow
  common generics;
- a full predictive-modelling toolchain: `logLik` / `AIC` / `BIC` /
  `perplexity` / `score_sequences` / `score_positions` / `model_fit`,
  five smoothing schemes (`floor`, `laplace`, `kneser_ney`,
  `witten_bell`, `jelinek_mercer`), k-fold cross-validated `tune_tree()`,
  permutation-tested two-tree `compare_trees()`, the multi-group
  `compare_groups()` (behavioral and usage axes, Benjamini-Hochberg FDR),
  and `bootstrap_pathways()` stability/informativeness with bootstrap CIs;
- sequence tools: `impute_sequences()` fills internal gaps from the
  fitted tree; `mine_contexts()` / `mine_sequences()` scan for contexts
  where a state is unusually likely or unlikely and for the best/worst-fit
  held-out sequences;
- flexible ingestion — a wide character matrix/data.frame, a list of
  character vectors, a long event log (via `actor` / `time` / `action` /
  `order` / `session` arguments, or the standalone `prepare_input()`), a
  `TraMineR` state-sequence exported with `as.matrix()` (the `%` void and
  `*` missing codes are treated as gaps), or a compatible `tna` /
  `Nestimate` transition/network object that carries its sequences — plus
  per-sequence weights for social-science workflows;
- five plot views — `horizontal` (default), `dendrogram`, `icicle`, an
  `interactive` HTML-widget tree, and the forward `plot_trajectories()` —
  with count-sized nodes and flow-sized edges.

## Installation

```r
# From CRAN
install.packages("transitiontrees")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("mohsaqr/transitiontrees")
```

## Usage

Every call below runs on the bundled datasets.

```r
library(transitiontrees)

# --- Fit ------------------------------------------------------------
# A wide character matrix / data.frame (one sequence per row), or a list
# of character vectors. `min_count` is the minimum occurrences a context
# needs to get its own node.
tree <- context_tree(trajectories, max_depth = 4, min_count = 5,
                     smoothing = "floor")

# Smoothing as an explicit list (hyperparameters), or a different scheme:
context_tree(trajectories, max_depth = 4, smoothing = list("floor", ymin = 0.001))
context_tree(trajectories, max_depth = 4, smoothing = "kneser_ney")

# From a long event log directly (no manual reshaping):
context_tree(group_regulation_long, actor = "Actor", time = "Time",
             action = "Action", max_depth = 3)
# or build the wide frame yourself:
wide <- prepare_input(group_regulation_long, actor = "Actor", time = "Time",
                      action = "Action")

# --- Inspect --------------------------------------------------------
print(tree)
summary(tree)
model_fit(tree)                          # logLik, df, nobs, AIC, BIC, perplexity

# --- Prune ----------------------------------------------------------
# Collapse contexts whose extra depth is not a significant gain over
# their parent (likelihood-ratio G^2 test at alpha = 0.05 per context).
pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)

# --- The pathway-centric API ---------------------------------------
tree_pathways(pruned)                       # all pathways, sorted by count
common_pathways(pruned,    top = 8)         # top by frequency
divergent_pathways(pruned, top = 6)         # top by KL from shorter history
divergent_pathways(pruned, flips_only = TRUE)
sharp_pathways(pruned,     top = 5)         # most deterministic continuations
tree_dependence(pruned)                     # per-context entropy/KL diagnostic

# --- Predict & evaluate --------------------------------------------
predict(pruned, newdata = list(c("Active", "Average"), c("Disengaged")))
predict(pruned, c("Active", "Active"), type = "class")   # modal next state
logLik(pruned); AIC(pruned); BIC(pruned)
perplexity(pruned, newdata = engagement)
score_sequences(pruned, newdata = engagement)

# Generate sequences from the fitted tree
generate_sequences(pruned, n = 100, length = 20)
simulate(pruned, nsim = 100, length = 20, seed = 1)   # S3 alias

# --- Cross-validated tuning ----------------------------------------
tg <- tune_tree(trajectories, max_depth = 1:4, min_count = c(3, 5),
                smoothing = c("floor", "kneser_ney"),
                prune = c(FALSE, TRUE), folds = 5)
tg                                          # ranked grid + the best configuration
plot(tg)

# --- Group fits and comparisons ------------------------------------
grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", group = "Achiever", max_depth = 2)
compare_trees(grp, iter = 200)              # two-tree permutation test
compare_groups(grp, iter = 199)             # per-context behavioral + usage tests

# --- Bootstrap pathway reliability ---------------------------------
boot <- bootstrap_pathways(pruned, iter = 1000, stat = "count")
summary(boot)        # tidy per-pathway stability + informativeness table
plot(boot)

# --- Tree introspection --------------------------------------------
query_pathway(pruned, "Active -> Average")   # full predicted distribution
pathway_exists(pruned, c("Active", "Average"))
subtree(pruned, "Active")                    # restrict to descendants of "Active"

# --- Impute & mine --------------------------------------------------
gappy <- list(c("Active", NA, "Average"), c("Average", "Average", NA, "Active"))
impute_sequences(pruned, gappy)                              # fill internal gaps
mine_contexts(pruned, state = "Disengaged", min_prob = 0.5)  # where Disengaged is likely
mine_sequences(pruned, engagement, which = "surprising")     # worst-fit sequences

# --- Visualisation --------------------------------------------------
plot(pruned)                                # horizontal tree (default)
plot(pruned, style = "dendrogram")          # pure-ggplot dendrogram
plot(pruned, style = "icicle")              # partition (sunburst) diagram
plot(pruned, style = "interactive")         # interactive HTML widget (visNetwork)
plot_pathways(pruned)                       # next-move probability heatmap
plot_divergence(pruned)                     # per-context KL lollipop
plot_distributions(pruned)                  # per-context next-state bars
plot_pruning(tree, "Active -> Active -> Average")   # suffix-chain memory view
plot_trajectories(pruned, measure = "frequency")    # forward trajectory tree
plot_predictive(pruned, engagement)         # held-out confidence diagnostics
plot_difference(grp, depth = 1)             # where two groups differ, per context
```

## Bundled datasets

Four datasets ship with the package for examples and tests:

| Dataset | Shape | Notes |
|---|---|---|
| `trajectories` | 138 × 15 wide character matrix, trailing `NA` = end | engagement states over time |
| `group_regulation_long` | 27,533-row long event log, POSIXct time | collaborative regulation moves, `Achiever` group label |
| `ai_long` | 8,551-row long event log, Unix time + session id | AI-collaboration moves |
| `engagement` | 1000 × 25 wide data.frame, `NA` = dropout | engagement-state sequences |

## Pathway notation

Pathways are reported in arrow notation (`A -> B -> C`). The leftmost
state is the *oldest*; the next-state prediction is conditional on the
trajectory ending at the rightmost state. The root context (the marginal
next-state distribution) is shown as `(start)`.

## Bootstrap interpretation

`bootstrap_pathways()` reports stability and informativeness separately.
`p_stability` is the bootstrap-estimated probability that a pathway
statistic falls outside the chosen consistency band under sequence-level
resampling. Small values mean the pathway rarely fails the reproducibility
criterion. `informative_rate` is the fraction of resamples where the
pathway's empirical G² statistic exceeds the chi-square reference cutoff
against its suffix-parent context.

Read the two flags together:

| stable | informative | Interpretation |
|---|---|---|
| TRUE | TRUE | reproducible and predictively distinctive pathway |
| TRUE | FALSE | reproducible pathway, but not predictively distinctive |
| FALSE | TRUE | sharp/divergent pathway carried by an unstable subset of sequences |
| FALSE | FALSE | weak or sample-fragile pathway |

## References

Begleiter, R., El-Yaniv, R., Yona, G. (2004). On prediction using
variable-order Markov models. *Journal of Artificial Intelligence
Research*, **22**, 385–421.

Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia: learning
probabilistic automata with variable memory length. *Machine Learning*,
**25**, 117–149.

Willems, F.M.J., Shtarkov, Y.M., Tjalkens, T.J. (1995). The context-tree
weighting method: basic properties. *IEEE Transactions on Information
Theory*, **41**, 653–664.

## License

GPL-3.
