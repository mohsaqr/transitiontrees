# transitiontrees v0.1.1 — Polish Pass

A function-by-function audit of the public surface for *dead-simple
calling, tidy outputs, and consistent conventions*. All changes are
backward-compatible at the **call** level except where noted.

## Inventory and verdicts

| Function | Audit | Change |
|----|----|----|
| [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md) | First-call OK (`context_tree(data)` works); 11 args feels heavy but each has a sensible default and most users will never touch the smoothing kwargs. | none |
| `prune()` | `prune(tree)` works (default G²). | none |
| [`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md) | Standard S3 `predict` shape; output matrix/named vector. | none |
| [`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md) | `n = 1` is a useless default — the function returns a 1×L matrix. | **default bumped to `n = 5`** |
| `pathways()` | Tidy data.frame; column names mostly OK but `length` was confusable (vs base R [`length()`](https://rdrr.io/r/base/length.html)) and the underlying concept is **depth**. | **column renamed `length` → `depth`** |
| [`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md) | Took a `length =` filter arg; same naming issue. | **arg renamed `length` → `depth`** |
| [`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md) / [`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md) | Wrap `pathways()`; inherit the rename automatically. | none |
| `path_dependence()` | Tidy data.frame, but column names diverged from `pathways()`: `context` vs `pathway`, `n` vs `count`. | **renamed `context` → `pathway`, `n` → `count`** |
| [`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md) / [`nobs.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/nobs.transitiontrees.md) | Standard `stats` generics; AIC / BIC work for free. | none |
| [`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md) | Took a `base =` arg that was *mathematically a no-op* (perplexity is base-invariant). Confusing. | **`base` arg removed** |
| [`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md) / [`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md) | Tidy data.frames. | none (already polished by `simplify` pass) |
| `smooth_transitiontrees()` | Many smoothing-method kwargs but each scheme genuinely needs its own. | none |
| `tune_transitiontrees()` | Returns `transitiontrees_tune` data.frame with `print` method, but no [`plot()`](https://rdrr.io/r/graphics/plot.default.html) — users couldn’t visualise their grid. | **added [`plot.transitiontrees_tune()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_tune.md)** |
| [`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md) / [`subtree()`](https://pak.dynasite.org/transitiontrees/reference/subtree.md) / [`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md) | Already simple; subtree adds `local_root` attribute. | none |
| `transitiontrees_distance()` | Bare scalar. | none |
| `compare_transitiontreess()` | Returns `transitiontrees_comparison` with `print` method, no [`plot()`](https://rdrr.io/r/graphics/plot.default.html). | **added [`plot.transitiontrees_comparison()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_comparison.md)** |
| [`plot.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md) | Two ggplot styles after the dependency cleanup. | none |
| [`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md) / [`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md) | Tidy lollipop / heatmap. | none |
| [`summary.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees.md) | Old summary table used `context` / `n` — same drift as `path_dependence`. | **table columns renamed to `pathway, depth, count, modal_next, prob_next`** |
| [`print.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/print.transitiontrees.md) | Old `(ctx) [n=23 p=0.5/0.3/0.2]` rendering was hard to scan; the per-node prob-vector slash list took the eye away from the modal next state — the thing users actually want at a glance. | **rewritten** — aligned columns, shows count + modal next + its probability, lists the smoothing scheme used |
| [`print.summary.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/print.summary.transitiontrees.md) | Same header re-styled. | banner format aligned with `print.transitiontrees` |
| **No `as.data.frame` method** | A user with a fitted tree shouldn’t have to call `pathways()` to get a flat node table — that’s `as.data.frame`’s job. | **added [`as.data.frame.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees.md)** returning `pathway, depth, count, modal_next, prob_next` |

## What “tidy output” means in transitiontrees

Every tidy data.frame returned by the package now uses the same column
vocabulary:

    pathway     character   arrow-notation pathway, e.g. "A -> B"; (root) for the root
    depth       integer     pathway length (depth in the tree)
    count       numeric     observed count of the pathway in training
    modal_next  character   the alphabet symbol with highest predicted probability
    prob_next   numeric     P(modal_next | pathway)
    KL          numeric     KL divergence vs. (k-1)-suffix; NA at root
    flips       logical     TRUE iff modal_next changes between this pathway and parent

These are the **canonical column names**. They appear in `pathways()`,
`path_dependence()`,
[`as.data.frame.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees.md),
and the
[`summary.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees.md)
table. They are stable.

## What “dead simple to call” means

Every public function has at most one required argument. All other
parameters have defaults that produce a useful result on a typical
dataset:

``` r

context_tree(data)                  # fits with sensible defaults
prune(tree)                         # G² with alpha = 0.05
pathways(tree)                      # all pathways, sorted by count
common_pathways(tree)               # top 10
divergent_pathways(tree)            # top 10 by KL
sharp_pathways(tree)                # top 10 by sharpness
path_dependence(tree)               # full diagnostic table
logLik(tree); AIC(tree); BIC(tree)  # standard model-comparison
perplexity(tree)                    # in-sample
generate_sequences(tree)            # 5 sequences of length 10
plot(tree)                          # dendrogram (ggplot)
as.data.frame(tree)                 # one-row-per-node tidy table
```

## What “every result has a plot()” means

| Result class | [`print()`](https://rdrr.io/r/base/print.html) | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) |
|----|----|----|
| `transitiontrees` | tree skeleton | dendrogram (default) / icicle (Suggests) |
| `summary.transitiontrees` | banner + node table | — |
| `transitiontrees_tune` | head + chosen line | **NEW**: perplexity surface, faceted by smoothing × prune, star = best |
| `transitiontrees_comparison` | observed + p + top divergent pathways | **NEW**: null-distribution histogram + observed line + p-value annotation |

## Backward compatibility

- **`pathways()$length` → `$depth`**: any code reading the `length`
  column directly will need to update.
  [`length()`](https://rdrr.io/r/base/length.html) (the base function)
  on the data.frame still works.
- **`common_pathways(length = k)` → `common_pathways(depth = k)`**: same
  rename of a kwarg.
- **`path_dependence()$context` / `$n` → `$pathway` / `$count`**.
- **`summary(tree)$table` columns**: `context, n, modal` →
  `pathway, count, modal_next` (plus new `prob_next`).
- **`perplexity(..., base = ...)`**: argument removed. Was a no-op.
- **`generate_sequences(n = 1)` → `n = 5`** by default.
- **[`print.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/print.transitiontrees.md)
  output** has changed visual format; existing scripts that grep its
  output may need updating.

## Tests

Tests for the renamed columns updated. Three new tests added
(`as.data.frame.transitiontrees`, `plot.transitiontrees_tune`,
`plot.transitiontrees_comparison`). **Suite: 273 tests, all pass.**
**PST 0.94.1 equivalence**: still machine precision (1.11e-16).
**`R CMD check`: 0 errors, 0 warnings, 0 notes.**
