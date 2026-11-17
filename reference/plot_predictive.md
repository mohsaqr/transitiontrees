# Plot how well a context tree predicts new sequences

`plot_predictive()` plots how much probability a context tree gives to
what happens next in sequences it has not seen, the measure of its value
as a predictive model. It scores every transition of held-out sequences
with
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
and shows where the model is confident and where it is surprised.

`type = "position"` plots the probability given to each observed next
state against its position in the sequence. `type = "ecdf"` plots the
empirical cumulative distribution of those probabilities, showing how
often the model gives high or low probability to what actually happened.
`type = "logloss"` plots the per-transition log-loss \\-\log_2
P(\mathrm{observed})\\ in bits against position, with the mean at each
position and the overall mean; 0 bits is a certain and correct
prediction.

For one number summarizing held-out prediction use
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md),
and to choose the settings that minimize it use
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

## Usage

``` r
plot_predictive(tree, newdata, type = c("position", "ecdf", "logloss"))
```

## Arguments

- tree:

  A `transitiontrees` object.

- newdata:

  Held-out sequences in any format accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- type:

  `"position"` (default), `"ecdf"` or `"logloss"`.

## Value

A ggplot object of the scored transitions.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object.

- `transitiontrees_empty_result`:

  (error) `newdata` has no transitions to score.

## See also

[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
for the scored transitions as a table;
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md);
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
held_out <- head(engagement, 40)
plot_predictive(tree, held_out)

plot_predictive(tree, held_out, type = "ecdf")

plot_predictive(tree, held_out, type = "logloss")
```
