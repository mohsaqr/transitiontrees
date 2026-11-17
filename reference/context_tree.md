# Fit a variable-order transition tree

`context_tree()` fits a variable-order transition tree to a set of
sequences. The tree, known in the literature as a prediction suffix tree
(Ron, Singer & Tishby, 1996) or a variable-length Markov chain (Bühlmann
& Wyner, 1999), describes how the next state of a process depends on the
states that precede it, while letting the amount of history that matters
differ from one history to another. Each node of the tree is a context,
or pathway: the most recent states, read from oldest to newest
(`"A -> B -> C"`). The root holds the marginal distribution of states;
its children condition on the last state, their children on the last two
states, and so on. The tree shows how much of the past is needed to
predict what comes next, and along which pathways a longer history
changes the prediction.

It counts every context of up to `max_depth` states that occurs at least
`min_count` times and estimates, for each, the conditional distribution
of the next state, smoothed so that unobserved continuations keep a
small probability. Sequences may be given as a wide table or matrix (one
row per sequence, one column per time point), a list of character
vectors, a long event log (one row per event; name the state column in
`action`), or a fitted network object from tna or Nestimate that carries
its sequences. With `group`, or with a grouped network object, one tree
is fitted per group over a shared set of states. With the default floor
smoothing, the fitted counts and probabilities match those of the PST
package (Gabadinho & Ritschard, 2016). It returns a `transitiontrees`
object, or a `transitiontrees_group` with one tree per group.

The tree is returned unpruned: every context that clears `min_count` is
kept. Use
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
to retain only the contexts whose prediction differs from that of their
shorter suffix,
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
to choose `max_depth`, `min_count` and the smoothing scheme by
cross-validation, and
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
or [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
read the fitted contexts as a table. Group differences are tested with
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
and the stability of individual pathways with
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).

## Usage

``` r
context_tree(
  data,
  max_depth = 5L,
  min_count = 5L,
  smoothing = "floor",
  alphabet = NULL,
  weights = NULL,
  group = NULL,
  block = NULL,
  actor = NULL,
  time = NULL,
  action = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900
)
```

## Arguments

- data:

  Sequences: a wide `data.frame` or character matrix (one row per
  sequence, one column per time point), a list of character vectors, a
  long event log (with `action`), or a tna or Nestimate model object
  that carries its sequence data. `NA`, `""` and the TraMineR codes
  `"%"` and `"*"` mark positions without an observation and are dropped.
  Numeric matrices are not accepted; convert them to character first.

- max_depth:

  Longest context, in states, that the tree may contain. Default 5.

- min_count:

  Minimum number of occurrences for a context to become a node. Default
  5.

- smoothing:

  Smoothing scheme for the next-state probabilities: a method name or a
  list of the method name and its hyperparameters. Methods are `"floor"`
  (default; `ymin = 0.001`, `rule = "interpolate"` or `"cap"`),
  `"laplace"` (`alpha = 1`), `"kneser_ney"` (`discount = 0.75`),
  `"witten_bell"` and `"jelinek_mercer"` (`lambda = 0.5`). See
  [`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md).

- alphabet:

  Character vector of states. `NULL` (default) uses the states observed
  in `data` (for a network object, its node labels). Supply it when new
  data may contain states unseen here.

- weights:

  Non-negative weight per sequence, one per row or list element of
  `data`. `NULL` (default) fits unweighted.

- group:

  Grouping for one tree per group: a vector with one value per sequence,
  the name of a column of `data` (long format) or of a network object's
  metadata. `NULL` (default) fits a single tree.

- block:

  Optional block (stratum) identifier in the same forms as `group`, kept
  with the fitted group and used by
  [`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
  for a stratified permutation test. Default `NULL`.

- actor, time, order, session, time_threshold:

  Long-format arguments passed to
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md):
  the column identifying whose events form a sequence, the timestamp
  column, an ordering column, a session column, and the gap in seconds
  that starts a new session (default 900). Each defaults to `NULL`.

- action:

  Name of the state column of a long event log. Giving it reshapes
  `data` with
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  before fitting. Default `NULL` (`data` is already in sequence form).

## Value

An object of class `transitiontrees`, the fitted tree. Read it through
its accessors:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) or
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
for one row per context,
[`state_frequencies()`](https://pak.dynasite.org/transitiontrees/reference/state_frequencies.md)
for the state counts,
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
for the log-likelihood, information criteria and perplexity, and
[`n_nodes()`](https://pak.dynasite.org/transitiontrees/reference/n_nodes.md)
for the number of contexts. With `group`, or a grouped input, an object
of class `transitiontrees_group`: one `transitiontrees` per group, in
group order, summarized by
[`summary()`](https://rdrr.io/r/base/summary.html) with one row per
group and flattened by
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to one
row per group and context.

## Conditions

`transitiontrees_bad_input` when `max_depth` or `min_count` is not a
single valid number, `data` has an unsupported type or no usable
sequence, a named column is absent, `weights`, `group` or `block` do not
have one entry per sequence, a network object carries no sequence data,
or the smoothing specification is invalid.

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

## See also

[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
and
[`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
to refine a fitted tree;
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
to choose its settings;
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
to read it;
[`plot.transitiontrees()`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md)
to draw it;
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
to compare groups.

## Examples

``` r
# Wide sequences: one learner per row, one time point per column
tree <- context_tree(trajectories, max_depth = 3)
tree
#> <transitiontrees>  38 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |-- Average   n=70     -> Active (0.53)
#> |   |   `-- Disengaged  n=15     -> Active (0.67)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   |   |-- Active    n=53     -> Average (0.55)
#> |   |   |-- Average   n=66     -> Active (0.58)
#> |   |   `-- Disengaged  n=12     -> Average (0.83)
#> |   `-- Disengaged  n=29     -> Active (0.52)
#> |       |-- Active    n=6      -> Active (0.83)
#> |       |-- Average   n=13     -> Active (0.54)
#> |       `-- Disengaged  n=9      -> Average (0.56)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   |-- Active    n=74     -> Average (0.55)
#> |   |   |-- Average   n=64     -> Active (0.45)
#> |   |   `-- Disengaged  n=10     -> Average (0.50)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.57)
#> |   |   |-- Average   n=248    -> Average (0.72)
#> |   |   `-- Disengaged  n=63     -> Average (0.65)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> |       |-- Active    n=7      -> Disengaged (0.57)
#> |       |-- Average   n=64     -> Disengaged (0.48) 
#> ... 12 more nodes (use as.data.frame(x) or summary(x))
as.data.frame(tree)
#>                                   pathway depth count likely_next
#> 1                                 (start)     0  1870     Average
#> 2                                 Average     1   751     Average
#> 3                                  Active     1   658      Active
#> 4                        Active -> Active     2   433      Active
#> 5                      Average -> Average     2   419     Average
#> 6                              Disengaged     1   325  Disengaged
#> 7              Active -> Active -> Active     3   316      Active
#> 8           Average -> Average -> Average     3   248     Average
#> 9                       Active -> Average     2   160     Average
#> 10                      Average -> Active     2   144      Active
#> 11               Disengaged -> Disengaged     2   139  Disengaged
#> 12                  Average -> Disengaged     2   134     Average
#> 13                  Disengaged -> Average     2   122     Average
#> 14 Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#> 15           Active -> Average -> Average     3    80     Average
#> 16            Active -> Active -> Average     3    74     Average
#> 17            Average -> Active -> Active     3    70      Active
#> 18           Average -> Average -> Active     3    66      Active
#> 19           Average -> Active -> Average     3    64      Active
#> 20       Average -> Disengaged -> Average     3    64  Disengaged
#> 21       Disengaged -> Average -> Average     3    63     Average
#> 22       Average -> Average -> Disengaged     3    62     Average
#> 23            Active -> Average -> Active     3    53     Average
#> 24    Disengaged -> Average -> Disengaged     3    46     Average
#> 25    Average -> Disengaged -> Disengaged     3    45  Disengaged
#> 26    Disengaged -> Disengaged -> Average     3    31     Average
#> 27                   Disengaged -> Active     2    29      Active
#> 28                   Active -> Disengaged     2    23      Active
#> 29        Active -> Average -> Disengaged     3    17      Active
#> 30         Disengaged -> Active -> Active     3    15      Active
#> 31        Average -> Disengaged -> Active     3    13      Active
#> 32        Disengaged -> Average -> Active     3    12     Average
#> 33        Disengaged -> Active -> Average     3    10     Average
#> 34     Disengaged -> Disengaged -> Active     3     9     Average
#> 35         Active -> Active -> Disengaged     3     9      Active
#> 36     Active -> Disengaged -> Disengaged     3     7  Disengaged
#> 37        Active -> Disengaged -> Average     3     7  Disengaged
#> 38         Active -> Disengaged -> Active     3     6      Active
#>    next_probability  divergence changes_prediction
#> 1         0.4347594          NA                 NA
#> 2         0.6098535 0.113562458              FALSE
#> 3         0.6975684 0.349487162               TRUE
#> 4         0.7852194 0.028601567              FALSE
#> 5         0.6778043 0.014667961              FALSE
#> 6         0.4830769 0.403065556               TRUE
#> 7         0.8354430 0.011491875              FALSE
#> 8         0.7177419 0.006242261              FALSE
#> 9         0.5187500 0.122825881              FALSE
#> 10        0.5000000 0.147275598              FALSE
#> 11        0.6762590 0.109773627              FALSE
#> 12        0.5000000 0.046118168               TRUE
#> 13        0.5163934 0.176186405              FALSE
#> 14        0.7037037 0.008106021              FALSE
#> 15        0.5750000 0.101876320              FALSE
#> 16        0.5540541 0.024836699              FALSE
#> 17        0.5285714 0.233006163              FALSE
#> 18        0.5757576 0.022257659              FALSE
#> 19        0.4531250 0.020767156               TRUE
#> 20        0.4843750 0.029405449               TRUE
#> 21        0.6507937 0.011442647              FALSE
#> 22        0.5161290 0.026852413              FALSE
#> 23        0.5471698 0.013380529               TRUE
#> 24        0.6086957 0.048844466              FALSE
#> 25        0.6000000 0.079870032              FALSE
#> 26        0.7085484 0.198404690              FALSE
#> 27        0.5172414 0.129904345              FALSE
#> 28        0.3913043 0.347809514               TRUE
#> 29        0.4117647 0.313295862               TRUE
#> 30        0.6666667 0.077077745              FALSE
#> 31        0.5384615 0.006062059              FALSE
#> 32        0.8318333 0.393334181               TRUE
#> 33        0.5000000 0.547770136              FALSE
#> 34        0.5555556 0.106026615               TRUE
#> 35        0.4444444 0.024571720              FALSE
#> 36        0.5714286 0.312801206              FALSE
#> 37        0.5707143 0.201375668               TRUE
#> 38        0.8318333 0.677361037              FALSE

# A long event log, split into sessions at gaps over 15 minutes
context_tree(group_regulation_long, actor = "Actor", time = "Time",
             action = "Action", max_depth = 2)
#> <transitiontrees>  87 nodes, depth <= 2, 9 states  [unpruned]
#>   alphabet : adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis
#>   fit on   : 2000 sequences, 27533 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=27533  -> consensus (0.25)
#> |-- adapt     n=509    -> consensus (0.47)
#> |   |-- consensus  n=27     -> cohesion (0.48)
#> |   |-- coregulate  n=28     -> consensus (0.50)
#> |   |-- discuss   n=259    -> consensus (0.47)
#> |   |-- emotion   n=6      -> consensus (0.83)
#> |   |-- monitor   n=16     -> consensus (0.50)
#> |   |-- plan      n=6      -> consensus (0.66)
#> |   `-- synthesis  n=140    -> consensus (0.48)
#> |-- cohesion  n=1695   -> consensus (0.50)
#> |   |-- adapt     n=130    -> consensus (0.54)
#> |   |-- cohesion  n=45     -> consensus (0.42)
#> |   |-- consensus  n=84     -> consensus (0.51)
#> |   |-- coregulate  n=61     -> consensus (0.50)
#> |   |-- discuss   n=173    -> consensus (0.54)
#> |   |-- emotion   n=843    -> consensus (0.47)
#> |   |-- monitor   n=71     -> consensus (0.53)
#> |   |-- plan      n=147    -> consensus (0.52)
#> |   `-- synthesis  n=20     -> consensus (0.55)
#> |-- consensus  n=6329   -> plan (0.40)
#> |   |-- adapt     n=222    -> plan (0.37)
#> |   |-- cohesion  n=795    -> plan (0.38)
#> |   |-- consensus  n=483    -> plan (0.39)
#> |   |-- coregulate  n=252    -> plan (0.40)
#> |   |-- discuss   n=1180   -> plan (0.41)
#> |   |-- emotion   n=830    -> plan (0.39) 
#> ... 61 more nodes (use as.data.frame(x) or summary(x))

# One tree per achievement group over a shared set of states
groups <- context_tree(group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2)
summary(groups)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC
#>   High  1000 13721      80         2  FALSE -22773 640 13721 46826 51643
#>    Low  1000 13812      82         2  FALSE -22968 656 13812 47248 52190
#>  perplexity
#>        5.26
#>        5.27
```
