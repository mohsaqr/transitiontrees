# A complete workflow

A variable-order Markov chain describes a process in which the amount of
history needed to predict the next state depends on the history itself.
The model is stored as a prediction suffix tree: every node is a
context, a run of recent states, with a distribution over the next
state, and a context is kept only where the data support it (Ron,
Singer, & Tishby, 1996; Bühlmann & Wyner, 1999). Fitting such a tree is
one step of an analysis that also has to decide how deep the tree may
grow, which of its contexts are reliable, how well it predicts data it
has not seen, and whether two groups of sequences differ.

This guide takes one event log through the whole analysis: reading the
data, fitting the tree, reading its pathways, smoothing and pruning,
choosing the depth by cross-validation, assessing fit, checking the
reliability of the pathways with the bootstrap, predicting next states,
and comparing two groups. Each step is treated briefly; the specialized
guides give the detail: [data
input](https://pak.dynasite.org/transitiontrees/articles/data-input.md),
[smoothing and
pruning](https://pak.dynasite.org/transitiontrees/articles/smoothing-and-pruning.md),
[model
selection](https://pak.dynasite.org/transitiontrees/articles/model-selection.md),
[pathway
reliability](https://pak.dynasite.org/transitiontrees/articles/pathway-reliability.md),
[prediction and
scoring](https://pak.dynasite.org/transitiontrees/articles/prediction-and-scoring.md),
[comparing
groups](https://pak.dynasite.org/transitiontrees/articles/comparing-groups.md)
and
[visualization](https://pak.dynasite.org/transitiontrees/articles/visualization.md).

``` r

library(transitiontrees)
```

## Data

`group_regulation_long` is an event log of collaborative regulation of
learning: each row is one regulation move made by a student working in a
small group. The data are synthetic, generated for the `tna` package
from its `group_regulation` data; they serve as a well-structured
example and are not evidence about students.

- **Actor:** `Actor`, the student; one sequence is built per student and
  session.
- **Time:** `Time`, the time stamp of the move (POSIXct), used to order
  the moves and to split sessions.
- **State:** `Action`, the regulation move, one of nine codes (`adapt`,
  `cohesion`, `consensus`, `coregulate`, `discuss`, `emotion`,
  `monitor`, `plan`, `synthesis`).
- **Group covariate:** `Achiever`, whether the student is a high or a
  low achiever; used in the group comparison.
- **Other columns:** `Group` (the collaboration group) and `Course`; not
  used here.

``` r

head(group_regulation_long)
#>   Actor Achiever Group Course                Time    Action
#> 1     1     High     1      A 2025-01-01 08:27:07  cohesion
#> 2     1     High     1      A 2025-01-01 08:35:20 consensus
#> 3     1     High     1      A 2025-01-01 08:42:18   discuss
#> 4     1     High     1      A 2025-01-01 08:50:00 synthesis
#> 5     1     High     1      A 2025-01-01 08:52:25     adapt
#> 6     1     High     1      A 2025-01-01 08:57:31 consensus
with(group_regulation_long,
     c(events = length(Actor), students = length(unique(Actor)),
       groups = length(unique(Group)), courses = length(unique(Course))))
#>   events students   groups  courses 
#>    27533     2000      200        3
```

The log holds 27,533 moves by 2,000 students in 200 groups and three
courses.

## Fitting the tree

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads the long log directly when the state column is named in `action`.
The moves of each `actor` are ordered by `time`, and a gap of more than
`time_threshold` seconds (default 900) starts a new session. The tree
below allows histories of up to three moves and gives a context its own
node when it occurs at least ten times.

``` r

tree <- context_tree(group_regulation_long,
                     actor = "Actor", time = "Time", action = "Action",
                     max_depth = 3, min_count = 10)
tree
#> <transitiontrees>  377 nodes, depth <= 3, 9 states  [unpruned]
#>   alphabet : adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis
#>   fit on   : 2000 sequences, 27533 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 10
#> (start)   n=27533  -> consensus (0.25)
#> |-- adapt     n=509    -> consensus (0.47)
#> |   |-- consensus  n=27     -> cohesion (0.48)
#> |   |-- coregulate  n=28     -> consensus (0.50)
#> |   |   `-- consensus  n=21     -> consensus (0.43)
#> |   |-- discuss   n=259    -> consensus (0.47)
#> |   |   |-- consensus  n=60     -> consensus (0.53)
#> |   |   |-- coregulate  n=37     -> consensus (0.43)
#> |   |   |-- discuss   n=48     -> consensus (0.52)
#> |   |   |-- emotion   n=14     -> consensus (0.50)
#> |   |   |-- monitor   n=25     -> consensus (0.40)
#> |   |   `-- plan      n=26     -> consensus (0.42)
#> |   |-- monitor   n=16     -> consensus (0.50)
#> |   `-- synthesis  n=140    -> consensus (0.48)
#> |       `-- discuss   n=107    -> consensus (0.45)
#> |-- cohesion  n=1695   -> consensus (0.50)
#> |   |-- adapt     n=130    -> consensus (0.54)
#> |   |   |-- consensus  n=13     -> consensus (0.61)
#> |   |   |-- discuss   n=65     -> consensus (0.55)
#> |   |   `-- synthesis  n=32     -> consensus (0.56)
#> |   |-- cohesion  n=45     -> consensus (0.42)
#> |   |   `-- emotion   n=24     -> consensus (0.46)
#> |   |-- consensus  n=84     -> consensus (0.51)
#> |   |   |-- cohesion  n=11     -> plan (0.45)
#> |   |   |-- discuss   n=13     -> consensus (0.61)
#> |   |   |-- emotion   n=11     -> consensus (0.45) 
#> ... 351 more nodes (use as.data.frame(x) or summary(x))
```

The banner reports 2,000 sequences, one per student, and 27,533 observed
moves. The unpruned tree has 377 contexts. The first line of the tree is
the root, `(start)`: the next-move distribution with no history, in
which `consensus` is the most likely move with probability 0.25. Every
deeper context is judged against the prediction of its parent.

The root distribution is the marginal frequency of each move.

``` r

state_frequencies(tree)
#>        state count proportion
#> 1  consensus  6797     0.2469
#> 2       plan  6623     0.2405
#> 3    discuss  4267     0.1550
#> 4    emotion  3075     0.1117
#> 5 coregulate  2133     0.0775
#> 6   cohesion  1839     0.0668
#> 7    monitor  1516     0.0551
#> 8  synthesis   729     0.0265
#> 9      adapt   554     0.0201
```

The moves are unevenly used: `consensus` and `plan` account for about a
quarter of the moves each, `adapt` and `synthesis` for 2% to 3%. A model
that always predicts `consensus` is correct a quarter of the time, so
the informative question is which histories change that default.

## Common and divergent pathways

The pathway tables give the contexts of the tree in one schema: the
pathway (oldest move on the left), its depth, its count, the most likely
next move and its probability, the Kullback-Leibler divergence in bits
of its next-move distribution from that of its parent, and whether the
most likely next move differs from the parent’s (`changes_prediction`).
[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md)
ranks the contexts by count.

``` r

common_pathways(tree, top = 8)
#>             pathway depth count likely_next next_probability divergence
#> 1           (start)     0 27533   consensus            0.247         NA
#> 2         consensus     1  6329        plan            0.396   0.334030
#> 3              plan     1  6157        plan            0.374   0.234143
#> 4           discuss     1  3951   consensus            0.321   0.555626
#> 5           emotion     1  2837    cohesion            0.325   0.555145
#> 6 consensus -> plan     2  2336        plan            0.375   0.000648
#> 7      plan -> plan     2  2108        plan            0.376   0.000432
#> 8        coregulate     1  1970     discuss            0.274   0.180891
#>   changes_prediction
#> 1                 NA
#> 2               TRUE
#> 3               TRUE
#> 4              FALSE
#> 5               TRUE
#> 6              FALSE
#> 7              FALSE
#> 8               TRUE
```

The most frequent histories are the single moves and the pairs formed by
`consensus` and `plan`. Their divergences are small at depth 2: after
`consensus -> plan` the next move is predicted almost exactly as after
`plan` alone (0.0006 bits).
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md)
ranks by divergence instead.

``` r

divergent_pathways(tree, top = 8)
#>                             pathway depth count likely_next next_probability divergence
#> 1 synthesis -> discuss -> consensus     3    10  coregulate            0.596      0.940
#> 2     consensus -> cohesion -> plan     3    12        plan            0.827      0.926
#> 3                         synthesis     1   652   consensus            0.463      0.892
#> 4    cohesion -> discuss -> emotion     3    10    cohesion            0.496      0.872
#> 5                             adapt     1   509   consensus            0.474      0.789
#> 6     monitor -> monitor -> discuss     3    12     discuss            0.249      0.783
#> 7 cohesion -> cohesion -> consensus     3    19        plan            0.366      0.756
#> 8  coregulate -> emotion -> monitor     3    13   consensus            0.382      0.752
#>   changes_prediction
#> 1               TRUE
#> 2              FALSE
#> 3              FALSE
#> 4               TRUE
#> 5              FALSE
#> 6               TRUE
#> 7              FALSE
#> 8               TRUE
```

Six of the eight most divergent contexts are three-move histories seen
10 to 19 times, close to the minimum of ten that the fit allows. After
`consensus -> cohesion -> plan`, for example, `plan` follows with
probability 0.83, against 0.36 after `cohesion -> plan`, but the
estimate rests on 12 occurrences. A large divergence on a small count is
the pattern that pruning and the bootstrap are designed to test: it may
reflect a real routine or a handful of sequences.

## Smoothing and pruning

Each node carries a smoothed distribution; the default `floor` smoother
leaves a fully observed distribution unchanged and gives every
unobserved next move a probability of 0.001. Pruning then removes the
contexts whose history does not change the prediction of the parent by
more than chance.
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md)
applies the four criteria to the same tree.

``` r

compare_pruning(tree)
#>   criterion n_nodes reduction_pct
#> 1        G2      23          93.9
#> 2        KL     368           2.4
#> 3       AIC      20          94.7
#> 4       BIC      10          97.3
```

The likelihood-ratio criterion `G2` keeps 23 of the 377 contexts; `AIC`
keeps 20 and `BIC` 10, while `KL`, at its default threshold of 0.005
nats, keeps almost all of them. Most of the grown depth is therefore not
supported. The tree pruned with `G2` at $`\alpha = 0.05`$ is used below.

``` r

pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)
pruned
#> <transitiontrees>  23 nodes, depth <= 3, 9 states  [pruned]
#>   alphabet : adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis
#>   fit on   : 2000 sequences, 27533 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 10
#>   pruned by: G2   alpha = 0.05
#> (start)   n=27533  -> consensus (0.25)
#> |-- adapt     n=509    -> consensus (0.47)
#> |-- cohesion  n=1695   -> consensus (0.50)
#> |   `-- cohesion  n=45     -> consensus (0.42)
#> |-- consensus  n=6329   -> plan (0.40)
#> |   |-- cohesion  n=795    -> plan (0.38)
#> |   |   `-- cohesion  n=19     -> plan (0.37)
#> |   `-- emotion   n=830    -> plan (0.39)
#> |       `-- emotion   n=58     -> plan (0.62)
#> |-- coregulate  n=1970   -> discuss (0.27)
#> |-- discuss   n=3951   -> consensus (0.32)
#> |   |-- adapt     n=29     -> adapt (0.24)
#> |   `-- coregulate  n=486    -> consensus (0.32)
#> |       `-- discuss   n=88     -> consensus (0.27)
#> |-- emotion   n=2837   -> cohesion (0.33)
#> |   |-- emotion   n=199    -> cohesion (0.35)
#> |   `-- plan      n=831    -> consensus (0.33)
#> |       `-- cohesion  n=33     -> cohesion (0.27)
#> |-- monitor   n=1433   -> discuss (0.38)
#> |-- plan      n=6157   -> plan (0.37)
#> |   `-- cohesion  n=221    -> plan (0.36)
#> |       `-- consensus  n=12     -> plan (0.83)
#> `-- synthesis  n=652    -> consensus (0.46)
```

The pruned tree keeps the root, the nine single moves and 13 longer
histories. [Smoothing and
pruning](https://pak.dynasite.org/transitiontrees/articles/smoothing-and-pruning.md)
gives the formulas of the five smoothers and four criteria.

## Choosing the depth by cross-validation

In-sample criteria reward deeper trees.
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
estimates the predictive quality of each configuration on data the tree
has not seen: it divides the sequences into `folds` groups, fits on all
but one, and computes the perplexity of the held-out group. The grid
below crosses three depths with pruned and unpruned trees.

``` r

tuning <- tune_tree(group_regulation_long,
                    actor = "Actor", time = "Time", action = "Action",
                    max_depth = 1:3, min_count = 10, prune = c(FALSE, TRUE),
                    folds = 5, seed = 1)
summary(tuning)
#> <transitiontrees_tune summary>  6 configurations, 0 with a failed fold
#> selected (minimum held-out perplexity):
#>  max_depth min_count                           smoothing prune perplexity n_nodes_avg
#>          1        10 floor(ymin=0.001, rule=interpolate) FALSE       5.46          10
as.data.frame(tuning, what = "grid")
#>   max_depth min_count                           smoothing prune logLik n_scored
#> 1         1        10 floor(ymin=0.001, rule=interpolate) FALSE -46738    27533
#> 2         1        10 floor(ymin=0.001, rule=interpolate)  TRUE -46738    27533
#> 3         2        10 floor(ymin=0.001, rule=interpolate)  TRUE -46763    27533
#> 4         3        10 floor(ymin=0.001, rule=interpolate)  TRUE -46884    27533
#> 5         2        10 floor(ymin=0.001, rule=interpolate) FALSE -47111    27533
#> 6         3        10 floor(ymin=0.001, rule=interpolate) FALSE -48782    27533
#>   perplexity n_nodes_avg folds_failed
#> 1       5.46        10.0            0
#> 2       5.46        10.0            0
#> 3       5.47        12.2            0
#> 4       5.49        23.4            0
#> 5       5.53        79.6            0
#> 6       5.88       344.4            0
```

The lowest held-out perplexity, 5.46, belongs to the first-order tree.
The pruned trees of depth 2 and 3 are within 0.03 of it, whereas the
unpruned depth-3 tree reaches 5.88: without pruning, the deeper contexts
fit noise that does not recur in new sequences. For predicting the next
move, one move of history is enough in these data; the longer contexts
that survive pruning add little predictive accuracy, although they may
still be of substantive interest. [Model
selection](https://pak.dynasite.org/transitiontrees/articles/model-selection.md)
covers the tuning grid, the information criteria and the stability of
the selection across fold assignments.

## Assessing fit

[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
returns the log likelihood, the number of free parameters `df` (the
number of contexts times the number of states minus one), the number of
predicted moves, AIC, BIC and perplexity in one row.

``` r

rbind(unpruned = model_fit(tree), pruned = model_fit(pruned))
#>          logLik   df  nobs   AIC    BIC perplexity
#> unpruned -45465 3016 27533 96962 121763       5.21
#> pruned   -46570  184 27533 93509  95022       5.43
```

The unpruned tree has the higher in-sample likelihood and the lower
in-sample perplexity (5.21 against 5.43), but it spends 3,016 parameters
to obtain them; both AIC and BIC favor the pruned tree with 184
parameters. The in-sample perplexity of the pruned tree is close to its
held-out value in the tuning grid (5.49), which indicates that it does
not overfit.

## Reliability of the pathways

A context retained by pruning has passed a test in this sample.
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
asks whether it would pass again: it resamples whole sequences with
replacement, refits the tree on each resample, and reports for each
pathway whether its count is stable and whether its likelihood ratio
statistic against its parent exceeds the critical value in most
resamples.

``` r

boot <- bootstrap_pathways(pruned, iter = 200, seed = 1)
boot
#> <transitiontrees_bootstrap>  200 resamples
#>   stability  : count in [0.50, 1.50] x observed, p < 0.05
#>   informative: G^2 > qchisq(0.95, df=k-1) = 15.51, threshold 0.80
#>   pathways   : 23 total, 22 stable, 14 informative, 13 both
#> 
#> top pathways (stable + informative first):
#>                           pathway depth count p_stability stability_rate stable
#>                         consensus     1  6329       0.005              1   TRUE
#>                              plan     1  6157       0.005              1   TRUE
#>                           discuss     1  3951       0.005              1   TRUE
#>                           emotion     1  2837       0.005              1   TRUE
#>                        coregulate     1  1970       0.005              1   TRUE
#>                          cohesion     1  1695       0.005              1   TRUE
#>                           monitor     1  1433       0.005              1   TRUE
#>                         synthesis     1   652       0.005              1   TRUE
#>                             adapt     1   509       0.005              1   TRUE
#>  discuss -> coregulate -> discuss     3    88       0.005              1   TRUE
#>  informative_rate informative mean_G2 ci_G2_lo ci_G2_hi
#>              1.00        TRUE  2945.0   2759.0   3139.1
#>              1.00        TRUE  2000.9   1858.1   2146.1
#>              1.00        TRUE  3047.2   2922.9   3217.9
#>              1.00        TRUE  2191.6   2020.7   2363.4
#>              1.00        TRUE   499.9    421.0    592.7
#>              1.00        TRUE   748.4    652.8    862.6
#>              1.00        TRUE   463.3    386.3    554.9
#>              1.00        TRUE   838.1    721.2    981.2
#>              1.00        TRUE   576.0    497.2    678.0
#>              0.94        TRUE    25.3     14.1     39.5
#> # ... 13 more pathways (use summary(x) for full table)
```

Of the 23 pathways of the pruned tree, 22 are stable and 14 informative;
13 are both: the nine single moves and four longer histories.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the bootstrap estimates of every pathway with percentile intervals; the
longer histories are:

``` r

subset(as.data.frame(boot, what = "summary"), depth > 1,
       select = c(pathway, count, stable, informative, informative_rate,
                  ci_G2_lo, ci_G2_hi))
#>                              pathway count stable informative informative_rate ci_G2_lo
#> 10  discuss -> coregulate -> discuss    88   TRUE        TRUE            0.940    14.11
#> 11   emotion -> emotion -> consensus    58   TRUE        TRUE            0.960    13.81
#> 12                  adapt -> discuss    29   TRUE        TRUE            0.875     9.90
#> 13              cohesion -> cohesion    45   TRUE        TRUE            0.835    10.43
#> 15                   plan -> emotion   831   TRUE       FALSE            0.130     2.76
#> 16              emotion -> consensus   830   TRUE       FALSE            0.615     6.56
#> 17             cohesion -> consensus   795   TRUE       FALSE            0.700     7.82
#> 18             coregulate -> discuss   486   TRUE       FALSE            0.390     6.13
#> 19                  cohesion -> plan   221   TRUE       FALSE            0.220     3.91
#> 20                emotion -> emotion   199   TRUE       FALSE            0.775     9.08
#> 21       cohesion -> plan -> emotion    33   TRUE       FALSE            0.765     7.70
#> 22     consensus -> cohesion -> plan    12   TRUE       FALSE            0.635     7.41
#> 23 cohesion -> cohesion -> consensus    19  FALSE        TRUE            0.870     9.76
#>    ci_G2_hi
#> 10     39.5
#> 11     45.3
#> 12     43.0
#> 13     47.9
#> 15     20.5
#> 16     31.8
#> 17     37.0
#> 18     29.1
#> 19     28.0
#> 20     45.4
#> 21     40.4
#> 22     30.4
#> 23     42.3
```

The four longer histories that are both stable and informative are
`discuss -> coregulate -> discuss`, `emotion -> emotion -> consensus`,
`adapt -> discuss` and `cohesion -> cohesion`; their likelihood ratio
statistic exceeds the critical value of 15.5 in 84% to 96% of the
resamples. These are the contexts on which a claim about longer memory
can rest. The divergent `consensus -> cohesion -> plan` is retained by
pruning and its count is stable, but its statistic exceeds the critical
value in only 64% of the resamples, so its departure from
`cohesion -> plan` is not reproducible.

[Pathway
reliability](https://pak.dynasite.org/transitiontrees/articles/pathway-reliability.md)
describes the two criteria and the distribution of the resampled
statistics.

## Predicting the next move

[`predict()`](https://rdrr.io/r/stats/predict.html) applies the
longest-suffix rule: the prediction for a history is the distribution of
the longest context in the tree that ends the history. With
`tidy = TRUE` it returns one row per history with the matched context,
the most likely next move and the full distribution.

``` r

histories <- list(c("plan", "consensus"),
                  c("discuss", "coregulate", "discuss"),
                  "monitor")
predict(pruned, newdata = histories, tidy = TRUE)
#>   id                          history                          context predicted
#> 1  1                plan -> consensus                        consensus      plan
#> 2  2 discuss -> coregulate -> discuss discuss -> coregulate -> discuss consensus
#> 3  3                          monitor                          monitor   discuss
#>   probability   adapt cohesion consensus coregulate discuss emotion monitor  plan
#> 1       0.396 0.00474   0.0149     0.082     0.1877   0.188  0.0727  0.0466 0.396
#> 2       0.271 0.12487   0.0010     0.271     0.1361   0.204  0.1474  0.0010 0.001
#> 3       0.375 0.01117   0.0558     0.159     0.0579   0.375  0.0907  0.0181 0.216
#>   synthesis
#> 1   0.00758
#> 2   0.11361
#> 3   0.01605
```

The first history is predicted from `consensus` alone, because the
pruned tree holds no longer context ending in it; the second matches the
full three-move context retained by pruning. [Prediction and
scoring](https://pak.dynasite.org/transitiontrees/articles/prediction-and-scoring.md)
covers scoring of whole sequences, imputation and the generation of new
sequences.

## Comparing high and low achievers

Naming the external grouping column in `group` fits one tree per group
over a shared alphabet. The comparison uses depth 2, so that most
contexts are observed often enough in both groups.

``` r

groups <- context_tree(group_regulation_long,
                       actor = "Actor", time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2, min_count = 10)
summary(groups)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC perplexity
#>   High  1000 13721      74         2  FALSE -22787 592 13721 46758 51214       5.26
#>    Low  1000 13812      77         2  FALSE -22984 616 13812 47199 51840       5.28
```

Each group contributes 1,000 sequences, with 13,721 and 13,812 moves;
the two trees have 74 and 77 contexts.
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
tests two kinds of difference by permuting the group labels of whole
sequences: whether, after the same history, the groups choose different
next moves (`behavioral`, a count-weighted Jensen-Shannon divergence),
and whether they reach the same histories at the same rates (`usage`, a
G-squared homogeneity statistic).

``` r

comparison <- compare_groups(groups, iter = 99, seed = 1)
summary(comparison)
#> Group comparison: High vs Low (80 contexts, 99 permutations)
#>   behavioral omnibus: JSD = 1772.142 bits, p = 0.010
#>   usage omnibus     : G2 = 1356.46, p = 0.010
#>   pathways with a behavioral difference (adjusted p < 0.05): 28
#>   pathways with a usage difference      (adjusted p < 0.05): 43
#>   pathways whose modal next state differs between groups: 26
```

Both omnibus statistics exceed every one of the 99 permuted values (p =
0.01, the smallest attainable with 99 permutations), so high and low
achievers differ both in what they do after a given history and in which
histories they reach. The contexts ranked by behavioral divergence, with
Benjamini-Hochberg adjusted p-values:

``` r

head(as.data.frame(comparison, what = "pathways"))
#>                   pathway depth count_total count_High count_Low modal_High modal_Low
#> 1        adapt -> monitor     2          16          4        12    discuss   discuss
#> 2    monitor -> synthesis     2          22         11        11  consensus     adapt
#> 3 synthesis -> coregulate     2          26          4        22    emotion consensus
#> 4        adapt -> discuss     2          29          4        25  consensus     adapt
#> 5       plan -> synthesis     2          11         10         1  consensus     adapt
#> 6    synthesis -> discuss     2          36          7        29  consensus synthesis
#>   flips jsd_bits jsd_p jsd_padj usage_g2 usage_p usage_padj
#> 1 FALSE    0.381  0.15    0.279 4.12e+00    0.05     0.0823
#> 2  TRUE    0.350  0.11    0.244 3.92e-04    0.95     0.9622
#> 3  TRUE    0.316  0.09    0.218 1.36e+01    0.01     0.0198
#> 4  TRUE    0.269  0.43    0.583 1.68e+01    0.01     0.0198
#> 5  TRUE    0.189  0.86    0.913 8.63e+00    0.02     0.0367
#> 6  TRUE    0.179  0.41    0.575 1.43e+01    0.01     0.0198
```

[Comparing
groups](https://pak.dynasite.org/transitiontrees/articles/comparing-groups.md)
explains the two statistics, the stratified permutation for nested data
and the plots of group differences.

## Reporting

A report of this analysis states the data and session rule, the fitted
and selected configurations, and the evidence for each claim, for
example:

> Regulation moves were modeled as a variable-order Markov chain with a
> prediction suffix tree (Ron, Singer, & Tishby, 1996). Moves were
> grouped into one sequence per student and session (gap \> 900 s),
> giving 2,000 sequences and 27,533 moves over nine codes. A tree of
> depth 3 (contexts observed at least ten times; floor smoothing,
> $`y_{\min} = 0.001`$) was pruned by likelihood-ratio tests at
> $`\alpha = 0.05`$, which retained 23 of 377 contexts. Five-fold
> cross-validation favored one move of history (held-out perplexity
> 5.46; pruned depth 3: 5.49). In 200 bootstrap resamples, 13 pathways
> were stable and informative, four of them longer than one move. High
> and low achievers differed in next-move choices given the same history
> and in the histories they reached (permutation tests with 99
> permutations, both p = 0.01).

## References

Bühlmann, P., & Wyner, A. J. (1999). Variable length Markov chains. *The
Annals of Statistics*, 27(2), 480–513.
<https://doi.org/10.1214/aos/1018031204>

Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
probabilistic automata with variable memory length. *Machine Learning*,
25(2–3), 117–149. <https://doi.org/10.1023/A:1026490906255>
