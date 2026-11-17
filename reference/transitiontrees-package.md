# transitiontrees: Transition Trajectories and Dynamics of Variable-Length Pathways or Sequences

`transitiontrees` models how categorical processes unfold: how one state
follows another in learning sessions, collaborative work, care pathways
or any other ordered record of events. Its model is the variable-order
transition tree, known as a prediction suffix tree (Ron, Singer &
Tishby, 1996) or variable-length Markov chain (Bühlmann & Wyner, 1999).
A first-order Markov chain predicts the next state from the current
state alone, and a chain of fixed order \\m\\ from the last \\m\\
states, with a number of parameters that grows exponentially in \\m\\. A
transition tree lets the length of the relevant history vary: it keeps a
longer context only where the data show that the extra history changes
what comes next.

Each node of the tree is a pathway, a run of recent states written from
oldest to newest (`"A -> B -> C"`), with the distribution of the state
that followed it. The package treats these pathways as the unit of
analysis. It reports which pathways occur, which of them change the
prediction relative to a shorter history, how stable they are under
resampling, and whether they differ between groups, and it uses the
fitted tree to predict, score, simulate and impute sequences. All
results are returned as data frames with one row per pathway, sequence
or group, and all plots are ggplot2 objects.

## Details

With the default floor smoothing, the fitted counts and probabilities
match those of the PST package (Gabadinho & Ritschard, 2016). The
smoothing schemes other than the floor, the pruning criteria,
cross-validated tuning, the bootstrap and the permutation tests are this
package's own.

## Data input

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
accepts wide sequence tables, lists of sequences, long event logs
(reshaped by
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)),
and transition network models from tna (Saqr et al., 2025) and Nestimate
that carry their sequences. The bundled datasets are
[`trajectories`](https://pak.dynasite.org/transitiontrees/reference/trajectories.md),
[`engagement`](https://pak.dynasite.org/transitiontrees/reference/engagement.md),
[`group_regulation_long`](https://pak.dynasite.org/transitiontrees/reference/group_regulation_long.md)
and
[`ai_long`](https://pak.dynasite.org/transitiontrees/reference/ai_long.md).

## Fitting and refining a tree

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
fits a tree, or one tree per group.
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
removes contexts that do not change the prediction,
[`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
re-estimates probabilities under another smoothing scheme, and
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md)
and
[`compare_smoothing()`](https://pak.dynasite.org/transitiontrees/reference/compare_smoothing.md)
tabulate the alternatives.
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
chooses depth, minimum count, smoothing and pruning by cross-validated
perplexity.

## Reading pathways

[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
(also [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html))
lists every pathway;
[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md)
and
[`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md)
rank them by frequency, divergence from the shorter history, and
certainty of the next state.
[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md)
adds entropies,
[`state_frequencies()`](https://pak.dynasite.org/transitiontrees/reference/state_frequencies.md)
gives the state counts, and
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md),
[`pathway_exists()`](https://pak.dynasite.org/transitiontrees/reference/pathway_exists.md),
[`subtree()`](https://pak.dynasite.org/transitiontrees/reference/subtree.md)
and
[`mine_contexts()`](https://pak.dynasite.org/transitiontrees/reference/mine_contexts.md)
look up individual pathways.

## Fit, prediction and scoring

[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
returns log-likelihood, AIC, BIC and perplexity, in-sample or on new
data (see also
[`logLik.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/logLik.transitiontrees.md),
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md)
and
[`n_nodes()`](https://pak.dynasite.org/transitiontrees/reference/n_nodes.md)).
[`predict.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/predict.transitiontrees.md)
predicts the next state,
[`score_sequences()`](https://pak.dynasite.org/transitiontrees/reference/score_sequences.md),
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
and
[`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md)
score observed sequences,
[`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
simulates new ones, and
[`impute_sequences()`](https://pak.dynasite.org/transitiontrees/reference/impute_sequences.md)
fills gaps.

## Reliability and comparison

[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
resamples sequences to assess the stability of each pathway.
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
and
[`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)
compare two trees, and
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
tests group differences pathway by pathway with permutation tests.

## Visualization

[`plot.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md)
draws the tree;
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
[`plot_divergence()`](https://pak.dynasite.org/transitiontrees/reference/plot_divergence.md)
and
[`plot_distributions()`](https://pak.dynasite.org/transitiontrees/reference/plot_distributions.md)
show pathways and their next-state distributions;
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md)
shows the pruning decision along one pathway;
[`plot_predictive()`](https://pak.dynasite.org/transitiontrees/reference/plot_predictive.md)
shows predictive performance on new data;
[`plot_trajectories()`](https://pak.dynasite.org/transitiontrees/reference/plot_trajectories.md)
draws the forward tree of sequence beginnings; and
[`plot_difference()`](https://pak.dynasite.org/transitiontrees/reference/plot_difference.md)
and
[`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md)
show group differences and bootstrap distributions.

## References

Bühlmann, P., & Wyner, A. J. (1999). Variable length Markov chains. *The
Annals of Statistics*, 27(2), 480-513.
[doi:10.1214/aos/1018031204](https://doi.org/10.1214/aos/1018031204)

Gabadinho, A., & Ritschard, G. (2016). Analyzing state sequences with
probabilistic suffix trees: The PST R package. *Journal of Statistical
Software*, 72(3).
[doi:10.18637/jss.v072.i03](https://doi.org/10.18637/jss.v072.i03)

Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
probabilistic automata with variable memory length. *Machine Learning*,
25(2-3), 117-149.
[doi:10.1023/A:1026490906255](https://doi.org/10.1023/A%3A1026490906255)

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK '25)* (pp. 351-361).
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

## See also

Useful links:

- <https://github.com/mohsaqr/transitiontrees>

- <https://pak.dynasite.org/transitiontrees/>

- Report bugs at <https://github.com/mohsaqr/transitiontrees/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

Authors:

- Mohammed Saqr <saqr@saqr.me>
  ([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

- Sonsoles López-Pernas <sonsoles.lopez@uef.fi>
  ([ORCID](https://orcid.org/0000-0002-9621-1392))
