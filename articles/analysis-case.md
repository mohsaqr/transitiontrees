# 2. A complete analysis case: collaborative-regulation sequences

This vignette runs one dataset all the way through and **reads the
numbers at each step** – not just *what* to call but *what the output
means* and *what not to over-read*. The data are the bundled
`group_regulation_long` event log: students’ collaborative
regulation-of-learning actions (`plan`, `monitor`, `consensus`,
`discuss`, …), one row per action, with a `High` / `Low` achievement
label per student.

The arc that emerges, stated up front so the sections connect:
regulation talk has **short memory** – the immediately preceding action
carries most of the predictive signal – a handful of two-action routines
reproducibly add to it, and **high and low achievers regulate
differently**, which the permutation test confirms.

## 1. The data

``` r

data(group_regulation_long)
nrow(group_regulation_long)
#> [1] 27533
head(group_regulation_long)
#>   Actor Achiever Group Course                Time    Action
#> 1     1     High     1      A 2025-01-01 08:27:07  cohesion
#> 2     1     High     1      A 2025-01-01 08:35:20 consensus
#> 3     1     High     1      A 2025-01-01 08:42:18   discuss
#> 4     1     High     1      A 2025-01-01 08:50:00 synthesis
#> 5     1     High     1      A 2025-01-01 08:52:25     adapt
#> 6     1     High     1      A 2025-01-01 08:57:31 consensus
sort(table(group_regulation_long$Action), decreasing = TRUE)
#> 
#>  consensus       plan    discuss    emotion coregulate   cohesion    monitor 
#>       6797       6623       4267       3075       2133       1839       1516 
#>  synthesis      adapt 
#>        729        554
```

The nine actions are very unevenly used: `consensus` and `plan`
dominate, `adapt` and `synthesis` are rare. That imbalance is the most
important fact about the corpus and it echoes through every result – a
model that just guesses `consensus` will look deceptively good, so the
interesting question is never “what is the modal next action” but
“*which histories overturn that default*”.

## 2. Fit

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads the long log directly: name the unit (`actor`), the clock
(`time`), and the state (`action`); it reshapes into one sequence per
session and fits. Sessions are split where the time gap is large.

``` r

tree <- context_tree(group_regulation_long,
                     actor = "Actor", time = "Time", action = "Action",
                     max_depth = 3L, min_count = 10L)
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

The banner reports the depth, the node count, the alphabet, and the
sequence/observation totals. The root line is the **null model**: the
next action given *no* history. Every deeper context has to beat that to
earn its place.

## 3. Inspect

``` r

summary(tree)
#> <transitiontrees summary>  377 nodes, depth <= 3, 9 states  [unpruned]
#> 
#>     pathway depth count likely_next next_probability divergence
#>     (start)     0 27533   consensus        0.2468674         NA
#>   consensus     1  6329        plan        0.3957971  0.3340303
#>        plan     1  6157        plan        0.3742082  0.2341431
#>     discuss     1  3951   consensus        0.3211845  0.5556255
#>     emotion     1  2837    cohesion        0.3253437  0.5551452
#>  coregulate     1  1970     discuss        0.2736041  0.1808906
#>    cohesion     1  1695   consensus        0.4979351  0.3149856
#>     monitor     1  1433     discuss        0.3754361  0.2283039
#>   synthesis     1   652   consensus        0.4630613  0.8917924
#>       adapt     1   509   consensus        0.4741100  0.7892874
#>  changes_prediction
#>                  NA
#>                TRUE
#>                TRUE
#>               FALSE
#>                TRUE
#>                TRUE
#>               FALSE
#>                TRUE
#>               FALSE
#>               FALSE
#> # ... 367 more rows (use as.data.frame(tree) for the full table)
model_fit(tree)
#>      logLik   df  nobs      AIC      BIC perplexity
#> 1 -45464.76 3016 27533 96961.51 121762.5   5.213661
```

Perplexity is the readable scalar: the effective number of equally
likely next actions. The uniform baseline is 9 (nine actions, no
knowledge); the fitted tree’s 5.21 says recent history collapses nine
possibilities to about 5.2. Real structure – but this is *in-sample* and
the tree is over-grown, so read it as an optimistic bound. Sections 6
and 7 give the honest figure.

## 4. The pathway tables

Three named verbs each fix a useful sort over the one canonical schema.

``` r

common_pathways(tree, top = 8)      # the highways
#>             pathway depth count likely_next next_probability   divergence
#> 1           (start)     0 27533   consensus        0.2468674           NA
#> 2         consensus     1  6329        plan        0.3957971 0.3340302623
#> 3              plan     1  6157        plan        0.3742082 0.2341431253
#> 4           discuss     1  3951   consensus        0.3211845 0.5556255176
#> 5           emotion     1  2837    cohesion        0.3253437 0.5551452430
#> 6 consensus -> plan     2  2336        plan        0.3754281 0.0006484903
#> 7      plan -> plan     2  2108        plan        0.3757116 0.0004321472
#> 8        coregulate     1  1970     discuss        0.2736041 0.1808905912
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

``` r

divergent_pathways(tree, top = 8)   # where adding history changes the prediction most
#>                             pathway depth count likely_next next_probability
#> 1 synthesis -> discuss -> consensus     3    10  coregulate        0.5956000
#> 2     consensus -> cohesion -> plan     3    12        plan        0.8268333
#> 3                         synthesis     1   652   consensus        0.4630613
#> 4    cohesion -> discuss -> emotion     3    10    cohesion        0.4965000
#> 5                             adapt     1   509   consensus        0.4741100
#> 6     monitor -> monitor -> discuss     3    12     discuss        0.2487500
#> 7 cohesion -> cohesion -> consensus     3    19        plan        0.3661053
#> 8  coregulate -> emotion -> monitor     3    13   consensus        0.3821538
#>   divergence changes_prediction
#> 1  0.9397936               TRUE
#> 2  0.9259734              FALSE
#> 3  0.8917924              FALSE
#> 4  0.8716723               TRUE
#> 5  0.7892874              FALSE
#> 6  0.7829439               TRUE
#> 7  0.7560378              FALSE
#> 8  0.7523873               TRUE
```

``` r

sharp_pathways(tree, top = 8)       # the most peaked next-action predictions
#>                             pathway depth count likely_next next_probability
#> 1     consensus -> cohesion -> plan     3    12        plan        0.8268333
#> 2 discuss -> coregulate -> cohesion     3    11   consensus        0.7217273
#> 3  coregulate -> coregulate -> plan     3    14        plan        0.7088571
#> 4   consensus -> adapt -> consensus     3    12        plan        0.6616667
#> 5  cohesion -> discuss -> synthesis     3    12   consensus        0.6616667
#> 6    emotion -> discuss -> cohesion     3    14   consensus        0.6380714
#> 7           discuss -> plan -> plan     3    11   consensus        0.6316364
#> 8   emotion -> emotion -> consensus     3    58        plan        0.6161034
#>   divergence changes_prediction
#> 1  0.9259734              FALSE
#> 2  0.5087962              FALSE
#> 3  0.6616111              FALSE
#> 4  0.3790714              FALSE
#> 5  0.3634509              FALSE
#> 6  0.2248026              FALSE
#> 7  0.5802868               TRUE
#> 8  0.2405929              FALSE
```

Read the divergent table in two layers. The very top rows can have large
`divergence` on a small `count` – a short history seen just over the
`min_count` floor that happened to resolve one way. Those are
small-sample mirages; the bootstrap in section 8 exists to disarm them.
The rows that *also* carry a large `count` are the well-supported
redirections worth quoting.

The sharp table teaches the same caution from the probability side: a
`next_probability` near 1 on a low `count` is a near-empty cell after
smoothing, not a law of behaviour. Sharpness **with** support is a rule;
sharpness without it is noise.

## 5. Per-context diagnostics

[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md)
is the information-theoretic decomposition the KL pruning rule
thresholds: per context, how many **bits** of next-action uncertainty
the extra history removes (`entropy_drop`), and whether it flips the
modal prediction.

``` r

tree_dependence(tree, sort_by = "entropy_drop", top = 8)
#>                                pathway depth count divergence  entropy
#> 1        consensus -> cohesion -> plan     3    12  0.9259734 0.885185
#> 2       cohesion -> discuss -> discuss     3    13  0.6501997 1.599492
#> 3    synthesis -> discuss -> consensus     3    10  0.9397936 1.358019
#> 4 coregulate -> synthesis -> consensus     3    14  0.3945507 1.440219
#> 5     coregulate -> coregulate -> plan     3    14  0.6616111 1.347984
#> 6    discuss -> coregulate -> cohesion     3    11  0.5087962 1.160730
#> 7        monitor -> discuss -> monitor     3    13  0.4218380 1.520958
#> 8          monitor -> cohesion -> plan     3    10  0.4461585 1.432416
#>   entropy_before entropy_drop likely_next likely_before changes_prediction
#> 1       2.289429    1.4042437        plan          plan              FALSE
#> 2       2.683111    1.0836194   consensus     consensus              FALSE
#> 3       2.344887    0.9868676  coregulate          plan               TRUE
#> 4       2.383187    0.9429680        plan          plan              FALSE
#> 5       2.276847    0.9288633        plan          plan              FALSE
#> 6       2.067552    0.9068226   consensus     consensus              FALSE
#> 7       2.401484    0.8805256     discuss       discuss              FALSE
#> 8       2.289429    0.8570129   consensus          plan               TRUE
```

A large `entropy_drop` with `changes_prediction = TRUE` is the most
valuable kind of context: it both sharpens *and* redirects. Watch for
**negative** `entropy_drop` – the longer history left the next action
*more* uncertain than its parent; that is the textbook signature of a
context pruning should remove.

## 6. Prune to the reliable tree

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

The pruned banner reports the surviving node count and criterion;
compare it to the unpruned `tree` from section 2. Each removed context
failed a likelihood-ratio G-squared test against its one-shorter parent:
the extra history did not explain enough added variation in the next
action to justify keeping it. That the tree collapses so far is itself a
finding – most of the grown depth was unsupported, and the durable
structure lives near the root.

## 7. Held-out predictive quality

The honest, out-of-sample estimate comes from cross-validation, which
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
runs at the sequence level over a `(max_depth, min_count, ...)` grid –
no hand-made train/test split. The in-sample perplexity is the
optimistic bound; the cross-validated winner is the figure to report.

``` r

model_fit(pruned)                                  # in-sample (optimistic)
#>      logLik  df  nobs      AIC      BIC perplexity
#> 1 -46570.36 184 27533 93508.72 95021.78   5.427279

tg <- tune_tree(group_regulation_long,
               actor = "Actor", time = "Time", action = "Action",
               max_depth = 1L:3L, min_count = 10L, folds = 5L, seed = 1L)
tg                                                 # ranked grid; "best" = the CV winner
#> <transitiontrees_tune>  6 configurations
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1   10 floor(ymin=0.001, rule=interpolate) FALSE -46738.34    27533
#>          1   10 floor(ymin=0.001, rule=interpolate)  TRUE -46738.34    27533
#>          2   10 floor(ymin=0.001, rule=interpolate)  TRUE -46763.30    27533
#>          3   10 floor(ymin=0.001, rule=interpolate)  TRUE -46884.31    27533
#>          2   10 floor(ymin=0.001, rule=interpolate) FALSE -47110.74    27533
#>          3   10 floor(ymin=0.001, rule=interpolate) FALSE -48782.29    27533
#>  perplexity n_nodes_avg folds_failed
#>    5.460492        10.0            0
#>    5.460492        10.0            0
#>    5.465445        12.2            0
#>    5.489519        23.4            0
#>    5.534851        79.6            0
#>    5.881286       344.4            0
#> 
#> best (min perplexity):
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1   10 floor(ymin=0.001, rule=interpolate) FALSE -46738.34    27533
#>  perplexity n_nodes_avg folds_failed
#>    5.460492          10            0
```

A cross-validated perplexity close to the in-sample value is the
signature of a well-pruned model that generalises; a large gap would say
*prune harder*.

[`mine_sequences()`](https://pak.dynasite.org/transitiontrees/reference/mine_sequences.md)
then surfaces the sessions the fitted model predicts worst – the
atypical regulation trajectories worth a closer look – and
[`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)
the individual moves it is most blindsided by:

``` r

wide <- prepare_input(group_regulation_long,
                     actor = "Actor", time = "Time", action = "Action")
mine_sequences(pruned, newdata = wide, which = "surprising", n = 5L)
#>   sequence_id n_scored   log_lik perplexity
#> 1        1559        2 -8.349842   65.03470
#> 2         446        2 -6.908739   31.63833
#> 3        1823        3 -9.619186   24.68993
#> 4        1323        3 -9.542743   24.06875
#> 5        1671        3 -9.140954   21.05177
score_positions(pruned, newdata = wide, worst = 5L)
#>   sequence_id position matched_context observed predicted_prob   log_lik
#> 1          69       22            plan    adapt   0.0009745006 -6.933585
#> 2         235       17            plan    adapt   0.0009745006 -6.933585
#> 3         974       20            plan    adapt   0.0009745006 -6.933585
#> 4        1227        7            plan    adapt   0.0009745006 -6.933585
#> 5        1424        3            plan    adapt   0.0009745006 -6.933585
```

## 8. Bootstrap reliability

[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
asked “which contexts pass a criterion *in this dataset*?”. The
bootstrap asks the stricter question – “which pass *reproducibly* under
resampling?” – and reports two flags. **`stable`**: the count
reproduces. **`informative`**: the G-squared against the parent
reproducibly clears the chi-square bar. A claim worth making is
**both**.

``` r

boot <- bootstrap_pathways(pruned, iter = 100L, stat = "count", seed = 1L)
boot
#> <transitiontrees_bootstrap>  100 resamples
#>   stability  : count in [0.50, 1.50] x observed, p < 0.05
#>   informative: G^2 > qchisq(0.95, df=k-1) = 15.51, threshold 0.80
#>   pathways   : 23 total, 21 stable, 14 informative, 13 both
#> 
#> top pathways (stable + informative first):
#>                           pathway depth count p_stability stability_rate stable
#>                         consensus     1  6329        0.01              1   TRUE
#>                              plan     1  6157        0.01              1   TRUE
#>                           discuss     1  3951        0.01              1   TRUE
#>                           emotion     1  2837        0.01              1   TRUE
#>                        coregulate     1  1970        0.01              1   TRUE
#>                          cohesion     1  1695        0.01              1   TRUE
#>                           monitor     1  1433        0.01              1   TRUE
#>                         synthesis     1   652        0.01              1   TRUE
#>                             adapt     1   509        0.01              1   TRUE
#>  discuss -> coregulate -> discuss     3    88        0.01              1   TRUE
#>  informative_rate informative  mean_G2 ci_G2_lo ci_G2_hi
#>              1.00        TRUE 2953.021 2734.523 3152.769
#>              1.00        TRUE 1999.111 1875.164 2148.083
#>              1.00        TRUE 3052.763 2925.134 3219.423
#>              1.00        TRUE 2201.118 2025.104 2349.157
#>              1.00        TRUE  502.777  428.421  600.107
#>              1.00        TRUE  750.700  662.513  854.702
#>              1.00        TRUE  460.401  391.198  557.597
#>              1.00        TRUE  840.484  724.436  972.445
#>              1.00        TRUE  580.569  508.107  674.837
#>              0.96        TRUE   25.679   15.060   42.570
#> # ... 13 more pathways (use summary(x) for full table)
```

``` r

head(summary(boot), 10)
#>                             pathway depth count likely_next next_probability
#> 1                         consensus     1  6329        plan        0.3957971
#> 2                              plan     1  6157        plan        0.3742082
#> 3                           discuss     1  3951   consensus        0.3211845
#> 4                           emotion     1  2837    cohesion        0.3253437
#> 5                        coregulate     1  1970     discuss        0.2736041
#> 6                          cohesion     1  1695   consensus        0.4979351
#> 7                           monitor     1  1433     discuss        0.3754361
#> 8                         synthesis     1   652   consensus        0.4662577
#> 9                             adapt     1   509   consensus        0.4774067
#> 10 discuss -> coregulate -> discuss     3    88   consensus        0.2727273
#>    divergence changes_prediction        G2 p_stability stability_rate stable
#> 1   0.3340303               TRUE 2930.7338  0.00990099              1   TRUE
#> 2   0.2341431               TRUE 1998.5086  0.00990099              1   TRUE
#> 3   0.5556255              FALSE 3043.2993  0.00990099              1   TRUE
#> 4   0.5551452               TRUE 2183.3402  0.00990099              1   TRUE
#> 5   0.1808906               TRUE  494.0122  0.00990099              1   TRUE
#> 6   0.3149856              FALSE  740.1435  0.00990099              1   TRUE
#> 7   0.2283039               TRUE  453.5394  0.00990099              1   TRUE
#> 8   0.9091915              FALSE  821.7854  0.00990099              1   TRUE
#> 9   0.8132687              FALSE  573.8618  0.00990099              1   TRUE
#> 10  0.1663148              FALSE   20.2894  0.00990099              1   TRUE
#>    informative_rate informative flip_consistency mean_count  sd_count
#> 1              1.00        TRUE             0.92    6342.99 109.57655
#> 2              1.00        TRUE             0.92    6161.66 130.96821
#> 3              1.00        TRUE             0.92    3958.31  81.97335
#> 4              1.00        TRUE             0.67    2836.33  61.00390
#> 5              1.00        TRUE             0.99    1972.85  53.81062
#> 6              1.00        TRUE             0.92    1696.68  42.70970
#> 7              1.00        TRUE             1.00    1430.31  40.39144
#> 8              1.00        TRUE             0.92     658.20  26.31357
#> 9              1.00        TRUE             0.92     511.36  20.96390
#> 10             0.96        TRUE             0.80      87.01  10.51838
#>    ci_count_lo ci_count_hi mean_next_probability sd_next_probability
#> 1     6115.925    6536.925             0.3959500         0.006381042
#> 2     5857.650    6366.175             0.3732453         0.006934701
#> 3     3838.900    4157.050             0.3202636         0.006807873
#> 4     2735.600    2952.575             0.3301679         0.006276168
#> 5     1875.275    2073.575             0.2737128         0.010252119
#> 6     1611.425    1767.675             0.4979885         0.011727875
#> 7     1356.850    1511.625             0.3754592         0.013158816
#> 8      609.850     716.000             0.4691459         0.018529368
#> 9      466.425     544.150             0.4788970         0.021643330
#> 10      68.475     109.050             0.2821353         0.041049294
#>    ci_next_probability_lo ci_next_probability_hi mean_divergence sd_divergence
#> 1               0.3810029              0.4070170       0.3358281   0.010620980
#> 2               0.3598353              0.3869965       0.2340653   0.008670898
#> 3               0.3061619              0.3330049       0.5564292   0.014373664
#> 4               0.3198049              0.3451925       0.5598475   0.021406786
#> 5               0.2546615              0.2903713       0.1837969   0.015985660
#> 6               0.4750867              0.5160186       0.3191616   0.020813645
#> 7               0.3547976              0.4010721       0.2322337   0.021001349
#> 8               0.4384477              0.5138353       0.9209807   0.058113430
#> 9               0.4372385              0.5240070       0.8192246   0.051354931
#> 10              0.2117708              0.3637233       0.2138777   0.056185496
#>    ci_divergence_lo ci_divergence_hi    mean_G2      sd_G2   ci_G2_lo
#> 1         0.3148910        0.3529989 2953.02132 106.335090 2734.52287
#> 2         0.2178176        0.2533216 1999.11080  79.507861 1875.16412
#> 3         0.5348465        0.5860827 3052.76257  81.500119 2925.13359
#> 4         0.5151091        0.5996423 2201.11768  91.903624 2025.10387
#> 5         0.1588430        0.2133030  502.77714  47.114738  428.42086
#> 6         0.2866160        0.3609771  750.69969  52.467915  662.51281
#> 7         0.1978436        0.2723805  460.40086  43.010781  391.19790
#> 8         0.8041894        1.0322695  840.48398  64.700502  724.43592
#> 9         0.7347244        0.9244516  580.56882  41.164021  508.10688
#> 10        0.1270018        0.3397060   25.67886   7.081377   15.06023
#>      ci_G2_hi
#> 1  3152.76930
#> 2  2148.08313
#> 3  3219.42286
#> 4  2349.15734
#> 5   600.10734
#> 6   854.70155
#> 7   557.59747
#> 8   972.44484
#> 9   674.83712
#> 10   42.57021
```

[`summary()`](https://rdrr.io/r/base/summary.html) sorts the trustworthy
(stable *and* informative) pathways first, so the top rows are the
defensible set. The two flags screen different failure modes. `stable`
alone keeps high-count noise pathways; `informative` alone could surface
a low-count borderline pathway whose sample G-squared is high by chance.
Their conjunction is the defensible set.

``` r

plot(boot)
#> `height` was translated to `width`.
```

![](analysis-case_files/figure-html/bootstrap-plot-1.png)

In the forest plot each bar is a 95% bootstrap interval on G-squared;
the dashed line is the chi-square critical value. A bar entirely to the
right is reproducibly informative; a bar straddling the line is not safe
to claim.

## 9. Do high and low achievers regulate differently?

Fit **one tree per group** in a single call with `group =`, then test
where the groups diverge with a permutation null. The grouping variable
is an *external* student attribute (`Achiever`), not derived from the
actions themselves – otherwise the comparison would be circular.

``` r

grp <- context_tree(group_regulation_long,
                    actor = "Actor", time = "Time", action = "Action",
                    group = "Achiever", max_depth = 2L, min_count = 10L)
cmp <- compare_groups(grp, iter = 199L, seed = 1L)
cmp
#> <transitiontrees_group_comparison>  2 groups, 199 permutations
#>   groups   : High, Low
#>   contexts : 80 tested
#> 
#> omnibus (permutation):
#>        axis                 statistic    value p_value
#>  behavioral count-weighted JSD (bits) 1772.142   0.005
#>       usage                   sum G^2 1356.465   0.005
#> 
#> top pathways by behavioral divergence (JSD):
#>                  pathway depth count_total modal_High modal_Low flips jsd_bits
#>         adapt -> monitor     2          16    discuss   discuss FALSE    0.381
#>     monitor -> synthesis     2          22  consensus     adapt  TRUE    0.350
#>  synthesis -> coregulate     2          26    emotion consensus  TRUE    0.316
#>         adapt -> discuss     2          29  consensus     adapt  TRUE    0.269
#>        plan -> synthesis     2          11  consensus     adapt  TRUE    0.189
#>     synthesis -> discuss     2          36  consensus synthesis  TRUE    0.179
#>   consensus -> synthesis     2          38  consensus     adapt  TRUE    0.170
#>       monitor -> monitor     2          23    discuss   discuss FALSE    0.167
#>  coregulate -> synthesis     2          32  consensus consensus FALSE    0.157
#>         adapt -> emotion     2          57  consensus  cohesion  TRUE    0.146
#>  jsd_padj usage_padj
#>     0.293      0.077
#>     0.242      0.980
#>     0.293      0.011
#>     0.581      0.011
#>     0.912      0.027
#>     0.581      0.011
#>     0.531      0.980
#>     0.533      0.904
#>     0.581      0.370
#>     0.267      0.011
#> # ... 70 more (use as.data.frame(x) for the full table)
```

The printed omnibus table reports two axes. **behavioral** is the
count-weighted Jensen-Shannon divergence (bits) between the groups’
next-action distributions, summed over shared contexts – “given the same
history, do the groups do different things next?”. **usage** is the
summed G-squared homogeneity statistic – “do they reach the same
contexts at the same rates?”. Each `p_value` comes from permuting the
group labels. Below it, the pathway table ranks contexts by behavioral
divergence with Benjamini-Hochberg adjusted p-values (`jsd_padj`,
`usage_padj`); `as.data.frame(cmp)` returns the full table.

Both omnibus p-values sit at the permutation floor, `1 / (iter + 1)`: no
relabelling of the students produced a difference as large as the
observed one, so high and low achievers do regulate differently. The
printed rows also show that no single context carries the behavioral
difference on its own after FDR correction – the difference is spread
thinly over many histories, which is why the omnibus statistic, not any
one pathway, is the claim to report.

``` r

plot_difference(grp, depth = 1L)
```

![](analysis-case_files/figure-html/groups-plot-1.png)

The per-context residual map shows *where* the groups differ: red and
blue cells are the contexts a high achiever and a low achiever resolve
toward different next actions. `depth = 1L` restricts it to the
single-action contexts so the rows stay readable; drop it (or raise it)
to inspect deeper histories.

## Synthesis

Pulling the thread through every section:

1.  **The action alphabet is imbalanced** – frequency is a misleading
    lens and modal predictions are trivially `consensus`/`plan`.
2.  **Memory is short** – pruning collapses the tree to a small set of
    contexts, and held-out perplexity confirms the shallow model
    generalises.
3.  **The insight is in the divergent, well-counted contexts** – not the
    common ones, and not the spectacular low-count tail.
4.  **Only the stable-and-informative pathways are claimable** – the
    bootstrap is the trust filter between an eyeballed table and a
    finding.
5.  **High and low achievers regulate measurably differently** – the
    permutation test licenses the claim that the omnibus statistic is
    real, not a relabelling artefact.

Each claim is anchored to a function whose output you can re-run – the
whole point of a pathway-centric, testable model.
