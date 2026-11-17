# Plot the pathways whose history changes the prediction most

`plot_divergence()` ranks pathways by the Kullback-Leibler divergence,
in bits, between their next-state distribution and that of their shorter
suffix, the context without the earliest state. In a variable-order
model a longer history matters only where it changes what is expected
next. Each pathway is a lollipop whose point size is its count; points
are colored by whether the most likely next state differs from the
suffix's, and those pathways are annotated with the new most likely
state and its probability. The chart identifies the histories that carry
information a shorter memory would lose.

The pathways are those returned by
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md).
A large divergence on a small count can be a sampling accident; check
reproducibility with
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
and the pruning decision for one pathway with
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md).

## Usage

``` r
plot_divergence(tree, top = 15L, min_count = 5L, title = NULL, ...)
```

## Arguments

- tree:

  A `transitiontrees` object.

- top:

  Number of pathways shown. Default 15.

- min_count:

  Pathways with a lower count are dropped. Default 5.

- title:

  Plot title. Default `NULL`: a title giving the number of pathways.

- ...:

  Ignored.

## Value

A ggplot object: one lollipop per pathway, divergence (bits) on the
x-axis.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object.

- `transitiontrees_empty_result`:

  (error) no pathway reaches `min_count`.

## See also

[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md)
for the table;
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md),
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
plot_divergence(tree, top = 10)
```
