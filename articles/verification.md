# Verification against the PST package

Two questions can be asked of statistical software. Verification asks
whether the program computes what its method defines: the counts, the
estimates, the likelihood. Validation asks whether the method, as
computed, describes the process it is applied to. This article is about
the first question for transitiontrees. The package is an independent
implementation of the variable-order prediction suffix tree (Ron, Singer
and Tishby, 1996), also known as the variable length Markov chain
(Bühlmann and Wyner, 1999). Its code shares nothing with the reference
implementation of the model in R, the PST package (Gabadinho and
Ritschard, 2016), which was archived on CRAN after version 0.94.1.

The verification has two parts. The first runs in this article: on a
small bundled dataset, every quantity the tree estimates is computed
again by hand with a few lines of base R and compared with the package’s
output. The second is an external equivalence suite that fits the same
models with PST 0.94.1, the markovchain package and the tna package and
compares the results quantity by quantity. The suite lives outside the
package, in the source repository, so that neither `R CMD check` nor a
user needs the archived PST package.

## What the tree estimates

For a context $`s`$, a sequence of up to `max_depth` states read oldest
to newest, let $`N(s, a)`$ be the number of times state $`a`$ follows
$`s`$ in the data and $`N(s) = \sum_a N(s, a)`$. Every position of every
sequence is counted once for the root and once for each preceding run of
states up to `max_depth`, so the first state of a sequence is counted
only at the root and the second at the root and at depth one. A context
enters the tree when $`N(s) \ge`$`min_count`. Its maximum-likelihood
next-state distribution is

``` math
\hat P(a \mid s) = N(s, a) / N(s).
```

The default `floor` smoother changes only a distribution in which some
state has a zero count. With $`k`$ states and floor $`y`$ (`ymin`), it
becomes

``` math
\tilde P(a \mid s) = (1 - k\,y)\,\hat P(a \mid s) + y,
```

which is the rule used by PST. To predict a position, the tree uses the
longest context it holds among the states that precede it, and the
log-likelihood of a set of sequences is the sum, over all positions, of
the log-probability the tree gives to the observed state.

## Data

The hand computation uses the first ten learners of `trajectories`.

- **Sequence:** one row per learner.
- **Time:** columns `1` to `15`; trailing `NA`s mark the end of a
  shorter trajectory.
- **State:** `Active`, `Average` or `Disengaged`.

``` r

library(transitiontrees)
learners <- head(trajectories, 10)
learners
#>       1         2            3            4            5            6            7        
#>  [1,] "Active"  "Disengaged" "Disengaged" "Disengaged" "Active"     "Active"     "Active" 
#>  [2,] "Average" "Average"    "Average"    "Average"    "Average"    "Average"    "Average"
#>  [3,] "Average" "Active"     "Active"     "Active"     "Active"     "Active"     "Active" 
#>  [4,] "Active"  "Active"     "Active"     "Active"     "Active"     "Average"    "Average"
#>  [5,] "Active"  "Active"     "Active"     "Average"    "Active"     "Average"    "Active" 
#>  [6,] "Average" "Average"    "Disengaged" "Average"    "Disengaged" "Disengaged" "Average"
#>  [7,] "Active"  "Average"    "Average"    "Average"    "Average"    "Average"    "Average"
#>  [8,] "Average" "Disengaged" "Average"    "Disengaged" "Disengaged" "Disengaged" "Active" 
#>  [9,] "Active"  "Active"     "Active"     "Active"     "Active"     "Active"     "Average"
#> [10,] "Average" "Average"    "Average"    "Average"    "Average"    "Average"    "Active" 
#>       8         9         10           11        12           13        14        15       
#>  [1,] "Active"  "Average" "Active"     "Average" "Active"     NA        NA        NA       
#>  [2,] "Average" "Average" "Active"     "Average" "Average"    "Average" "Average" "Average"
#>  [3,] "Active"  "Active"  "Average"    "Active"  "Active"     "Active"  "Active"  "Active" 
#>  [4,] "Active"  "Active"  "Active"     "Average" "Active"     "Average" "Active"  "Average"
#>  [5,] "Active"  "Active"  "Disengaged" "Active"  "Active"     "Active"  "Average" "Average"
#>  [6,] "Average" "Average" "Disengaged" "Average" "Disengaged" "Average" "Average" "Average"
#>  [7,] "Average" "Average" "Disengaged" "Average" "Average"    "Average" "Average" "Average"
#>  [8,] "Active"  "Active"  "Active"     "Active"  "Active"     "Active"  "Active"  "Active" 
#>  [9,] "Active"  "Average" "Average"    "Active"  "Active"     "Active"  "Active"  "Active" 
#> [10,] "Average" "Average" "Disengaged" "Average" "Average"    "Average" "Average" "Average"
```

## Counting by hand

The hand computation needs two small functions written in base R, which
use nothing from the package. `count_by_hand()` lists, for every
position of every sequence, the root and each preceding run of up to
`max_depth` states, and tabulates the state that follows. It returns one
row per context with its total count and its count of each next state.
`contexts_of()` turns the pathway names of that table into histories
that can be passed to
[`predict()`](https://rdrr.io/r/stats/predict.html).

``` r

count_by_hand <- function(sequences, max_depth) {
  events <- lapply(seq_len(nrow(sequences)), function(i) {
    states <- unname(sequences[i, ])
    states <- states[!is.na(states)]
    per_depth <- lapply(seq.int(0L, max_depth), function(d) {
      positions <- seq_along(states)
      positions <- positions[positions > d]
      pathway <- vapply(positions, function(t) {
        if (d == 0L) "(start)" else
          paste(states[seq.int(t - d, t - 1L)], collapse = " -> ")
      }, character(1))
      data.frame(pathway = pathway, next_state = states[positions])
    })
    do.call(rbind, per_depth)
  })
  events <- do.call(rbind, events)
  counts <- table(events$pathway, events$next_state)
  data.frame(pathway = rownames(counts),
             count = as.integer(rowSums(counts)),
             as.data.frame.matrix(counts), row.names = NULL)
}

contexts_of <- function(hand) {
  lapply(hand$pathway, function(p) {
    if (identical(p, "(start)")) character(0) else
      strsplit(p, " -> ", fixed = TRUE)[[1L]]
  })
}
```

To count the contexts up to depth two, we call `count_by_hand()`.

``` r

hand <- count_by_hand(learners, max_depth = 2L)
hand
#>                     pathway count Active Average Disengaged
#> 1                   (start)   147     64      68         15
#> 2                    Active    60     43      15          2
#> 3          Active -> Active    40     32       7          1
#> 4         Active -> Average    14      8       6          0
#> 5      Active -> Disengaged     2      1       0          1
#> 6                   Average    62     13      41          8
#> 7         Average -> Active    12      5       7          0
#> 8        Average -> Average    36      4      28          4
#> 9     Average -> Disengaged     8      0       6          2
#> 10               Disengaged    15      3       7          5
#> 11     Disengaged -> Active     3      3       0          0
#> 12    Disengaged -> Average     7      0       4          3
#> 13 Disengaged -> Disengaged     5      2       1          2
```

The ten learners contribute 147 observed states, all counted at the
root. All nine pairs of states occur as contexts of depth two, so the
hand count has 13 contexts.

## Counts

To fit the tree that keeps every observed context, we call
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
with `min_count = 1L`; `smoothing = list("floor", ymin = 0)` switches
the smoother off so that the probabilities are the raw
maximum-likelihood estimates.

``` r

tree <- context_tree(learners, max_depth = 2L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
tree
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 10 sequences, 147 observations
#>   smoothing: floor(ymin=0, rule=interpolate)   min_count = 1
#> (start)   n=147    -> Average (0.46)
#> |-- Active    n=60     -> Active (0.72)
#> |   |-- Active    n=40     -> Active (0.80)
#> |   |-- Average   n=12     -> Average (0.58)
#> |   `-- Disengaged  n=3      -> Active (1.00)
#> |-- Average   n=62     -> Average (0.66)
#> |   |-- Active    n=14     -> Active (0.57)
#> |   |-- Average   n=36     -> Average (0.78)
#> |   `-- Disengaged  n=7      -> Average (0.57)
#> `-- Disengaged  n=15     -> Average (0.47)
#>     |-- Active    n=2      -> Active (0.50)
#>     |-- Average   n=8      -> Average (0.75)
#>     `-- Disengaged  n=5      -> Active (0.40)
```

To compare the counts, we merge the hand table with the tree’s table
from [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) by
pathway.

``` r

counts <- merge(subset(hand, select = c(pathway, count)),
                subset(as.data.frame(tree), select = c(pathway, count)),
                by = "pathway", suffixes = c("_hand", "_tree"))
counts
#>                     pathway count_hand count_tree
#> 1                   (start)        147        147
#> 2                    Active         60         60
#> 3          Active -> Active         40         40
#> 4         Active -> Average         14         14
#> 5      Active -> Disengaged          2          2
#> 6                   Average         62         62
#> 7         Average -> Active         12         12
#> 8        Average -> Average         36         36
#> 9     Average -> Disengaged          8          8
#> 10               Disengaged         15         15
#> 11     Disengaged -> Active          3          3
#> 12    Disengaged -> Average          7          7
#> 13 Disengaged -> Disengaged          5          5
c(by_hand = nrow(hand), in_tree = n_nodes(tree), matched = nrow(counts))
#> by_hand in_tree matched 
#>      13      13      13
with(counts, identical(count_hand, count_tree))
#> [1] TRUE
```

The tree holds the same 13 contexts as the hand count, and every count
is identical. The number of observations the tree reports,
[`nobs()`](https://rdrr.io/r/stats/nobs.html), equals the root count.

``` r

nobs(tree)
#> [1] 147
```

## Maximum-likelihood probabilities

To compute the maximum-likelihood distributions by hand, we divide each
next state’s count by the context’s count. To obtain the tree’s
distributions, we call
[`predict()`](https://rdrr.io/r/stats/predict.html) with each context as
a history and `tidy = TRUE`.

``` r

mle_by_hand <- function(hand) {
  states <- setdiff(names(hand), c("pathway", "count"))
  data.frame(pathway = hand$pathway, hand[states] / hand$count)
}
by_hand <- mle_by_hand(hand)
fitted <- predict(tree, newdata = contexts_of(hand), tidy = TRUE)
subset(fitted, select = c(history, context, Active, Average, Disengaged))
#>                     history                  context Active Average Disengaged
#> 1                   (start)                  (start) 0.4354  0.4626    0.10204
#> 2                    Active                   Active 0.7167  0.2500    0.03333
#> 3          Active -> Active         Active -> Active 0.8000  0.1750    0.02500
#> 4         Active -> Average        Active -> Average 0.5714  0.4286    0.00000
#> 5      Active -> Disengaged     Active -> Disengaged 0.5000  0.0000    0.50000
#> 6                   Average                  Average 0.2097  0.6613    0.12903
#> 7         Average -> Active        Average -> Active 0.4167  0.5833    0.00000
#> 8        Average -> Average       Average -> Average 0.1111  0.7778    0.11111
#> 9     Average -> Disengaged    Average -> Disengaged 0.0000  0.7500    0.25000
#> 10               Disengaged               Disengaged 0.2000  0.4667    0.33333
#> 11     Disengaged -> Active     Disengaged -> Active 1.0000  0.0000    0.00000
#> 12    Disengaged -> Average    Disengaged -> Average 0.0000  0.5714    0.42857
#> 13 Disengaged -> Disengaged Disengaged -> Disengaged 0.4000  0.2000    0.40000
identical(with(fitted, context), with(by_hand, pathway))
#> [1] TRUE
all.equal(subset(fitted, select = c(Active, Average, Disengaged)),
          subset(by_hand, select = c(Active, Average, Disengaged)),
          tolerance = 1e-12)
#> [1] TRUE
```

Every history is predicted from its own context, without back-off, and
the probabilities agree with the hand computation to within
$`10^{-12}`$.

## Floor smoothing

Six of the hand-counted distributions contain a zero, for example
`Active -> Disengaged` (one `Active`, one `Disengaged`, no `Average`).
To check the default smoother, we apply the floor rule by hand with
`ymin = 0.001` and compare it with a tree fitted with the default
smoothing.

``` r

floor_by_hand <- function(probabilities, ymin) {
  states <- setdiff(names(probabilities), "pathway")
  p <- as.matrix(probabilities[states])
  k <- length(states)
  has_zero <- rowSums(p == 0) > 0
  p[has_zero, ] <- (1 - k * ymin) * p[has_zero, ] + ymin
  data.frame(pathway = probabilities$pathway, p)
}
smoothed_tree <- context_tree(learners, max_depth = 2L, min_count = 1L)
smoothed <- predict(smoothed_tree, newdata = contexts_of(hand), tidy = TRUE)
all.equal(subset(smoothed, select = c(Active, Average, Disengaged)),
          subset(floor_by_hand(by_hand, ymin = 0.001),
                 select = c(Active, Average, Disengaged)),
          tolerance = 1e-12)
#> [1] TRUE
```

## The log-likelihood

To compute the log-likelihood by hand, we score every position with the
hand-computed distribution of its longest preceding context and sum the
logarithms.

``` r

loglik_by_hand <- function(sequences, max_depth, probabilities) {
  states_all <- setdiff(names(probabilities), "pathway")
  p <- as.matrix(probabilities[states_all])
  rownames(p) <- probabilities$pathway
  per_sequence <- vapply(seq_len(nrow(sequences)), function(i) {
    states <- unname(sequences[i, ])
    states <- states[!is.na(states)]
    sum(vapply(seq_along(states), function(t) {
      d <- min(t - 1L, max_depth)
      context <- if (d == 0L) "(start)" else
        paste(states[seq.int(t - d, t - 1L)], collapse = " -> ")
      log(p[context, states[t]])
    }, numeric(1)))
  }, numeric(1))
  sum(per_sequence)
}
loglik_by_hand(learners, max_depth = 2L, probabilities = by_hand)
#> [1] -99.95
logLik(tree)
#> 'log Lik.' -99.95 (df=26)
all.equal(loglik_by_hand(learners, max_depth = 2L, probabilities = by_hand),
          as.numeric(logLik(tree)), tolerance = 1e-12)
#> [1] TRUE
```

## The same check on the full data

The same comparison can be run on all the bundled sequences at several
depths. `verify_by_hand()` combines the steps above and returns one row
per dataset and depth: the number of contexts found by hand and in the
tree, and the largest absolute differences in counts, probabilities and
log-likelihood.

``` r

verify_by_hand <- function(sequences, max_depth, label) {
  hand <- count_by_hand(sequences, max_depth)
  tree <- context_tree(sequences, max_depth = max_depth, min_count = 1L,
                       smoothing = list("floor", ymin = 0))
  tree_counts <- as.data.frame(tree)
  matched <- merge(hand, tree_counts, by = "pathway",
                   suffixes = c("_hand", "_tree"))
  by_hand <- mle_by_hand(hand)
  states <- setdiff(names(by_hand), "pathway")
  fitted <- predict(tree, newdata = contexts_of(hand), tidy = TRUE)
  data.frame(
    data = label, max_depth = max_depth,
    contexts_hand = nrow(hand), contexts_tree = n_nodes(tree),
    max_count_error = max(abs(matched$count_hand - matched$count_tree)),
    max_prob_error = max(abs(as.matrix(fitted[states]) -
                               as.matrix(by_hand[states]))),
    loglik_error = abs(loglik_by_hand(sequences, max_depth, by_hand) -
                         as.numeric(logLik(tree))))
}
checks <- Map(verify_by_hand,
              sequences = list(trajectories, trajectories, trajectories,
                               as.matrix(engagement), as.matrix(engagement)),
              max_depth = c(1L, 2L, 3L, 1L, 2L),
              label = c("trajectories", "trajectories", "trajectories",
                        "engagement", "engagement"))
do.call(rbind, checks)
#>           data max_depth contexts_hand contexts_tree max_count_error max_prob_error loglik_error
#> 1 trajectories         1             4             4               0              0    0.000e+00
#> 2 trajectories         2            13            13               0              0    2.274e-13
#> 3 trajectories         3            40            40               0              0    0.000e+00
#> 4   engagement         1             4             4               0              0    0.000e+00
#> 5   engagement         2            13            13               0              0    1.819e-12
```

On the 136 non-empty trajectories and the 1,000 engagement sequences, at
every depth, the tree holds exactly the contexts found by hand, the
counts and the probabilities are identical, and the log-likelihoods
differ by at most $`2 \times 10^{-11}`$, the rounding error of summing
thousands of logarithms in a different order.
([`as.matrix()`](https://rdrr.io/r/base/matrix.html) gives the
engagement data frame the same matrix form as `trajectories`, which the
hand-count function expects.)

## The external equivalence suite

The hand computation verifies the estimator against its definition. The
external suite verifies it against implementations written by other
authors. It is kept in `local_testing_and_equivalence/` in the source
repository, outside the built package, and runs only when the
environment variable `TRANSITREES_EQUIV_TESTS` is `true`. Its runner
installs PST 0.94.1 from the CRAN archive, installs the current source
of transitiontrees, and runs three test files.

**PST.** Thirty configurations are drawn at random with a fixed seed.
Each simulates first-order Markov sequences from a random transition
matrix with 2 to 4 states, 50, 80 or 120 sequences of length 15, 20 or
25, and fits both packages with `max_depth` from 2 to 4, `min_count`
(`nmin` in PST) of 3, 5 or 8, and `ymin` of 0, 0.001 or 0.005. TraMineR
builds the sequence object that PST requires. Every configuration is
compared on:

| transitiontrees | PST 0.94.1 | Criterion |
|----|----|----|
| node counts, all depths ([`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)) | `pstree()` node counts | identical integers |
| node probabilities, all depths | `pstree()` node probabilities | absolute difference \< $`10^{-10}`$ |
| set of contexts | `nodenames()` | identical sets |
| [`nobs()`](https://rdrr.io/r/stats/nobs.html) | [`nobs()`](https://rdrr.io/r/stats/nobs.html) | identical |
| [`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md) | `query()` | absolute difference \< $`10^{-10}`$ |
| per-position probabilities ([`score_positions()`](https://pak.dynasite.org/transitiontrees/reference/score_positions.md)) | `predict(decomp = TRUE)` | absolute difference \< $`10^{-9}`$ |
| [`logLik()`](https://rdrr.io/r/stats/logLik.html) | [`logLik()`](https://rdrr.io/r/stats/logLik.html) | relative difference \< $`10^{-10}`$ |
| [`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md) distributions and counts | `query()` and node counts | \< $`10^{-10}`$; identical |
| `plot_predictive(type = "logloss")` | $`-\log_2`$ of `predict(decomp = TRUE)` | absolute difference \< $`10^{-9}`$ |
| `prune_tree(criterion = "G2")` | `prune(gain = "G2")` | measured agreement (below) |

**markovchain.** A tree of depth one without smoothing is a first-order
Markov chain. Its transition probabilities are compared with the
maximum-likelihood estimate of `markovchain::markovchainFit()` on three
simulated datasets with 2, 3 and 4 states (absolute difference \<
$`10^{-12}`$).

**tna.**
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md),
which turns a long event log into sequences, is compared with
[`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html)
over every combination of its arguments (actor, time, order, action and
the session time threshold) on three long-format datasets from the
Nestimate package, and an explicit session column is checked to split
sequences by actor and session.

### Results

The suites were last run on 2026-11-17 against transitiontrees 0.1.5 on
R 4.5.2 (Linux), with PST 0.94.1 installed from the CRAN archive,
markovchain from CRAN and tna 1.5.8.

| Suite | Tests (expectations) | Failures | Max abs count error | Max abs probability error |
|----|----|----|----|----|
| PST 0.94.1 | 10 (296) | 0 | 0 | 0 (log-loss: 4.4e-16) |
| markovchain | 1 (1) | 0 | not applicable | 0 |
| tna | 2 (422) | 0 | not applicable | not applicable |

The PST suite fits 30 random configurations (2 to 4 states, depth 2 to
4, `min_count` of 3, 5 or 8 and `ymin` of 0, 0.001 or 0.005). Every
shared context has the same counts and probabilities, and the
per-position predictive probabilities agree exactly; the only non-zero
difference is a log-loss difference of 4.4e-16, at the level of
floating-point rounding. For the tna comparison the contract is that
both functions produce the same set of sequences for each of the 420
argument combinations. Since tna 1.3.1, sessions built from a time
column are ordered numerically by actor, while
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
orders them as text, so the rows agree in order for 153 of the 420
combinations and as a set for all of them. Row order does not change a
fitted tree.

### What is not at parity

Some functions implement the package’s own design and have no
counterpart in PST, so they are verified by the package’s unit tests
only.

- **Pruning.** `prune_tree(criterion = "G2")` and PST’s
  `prune(gain = "G2")` use the same likelihood-ratio statistic but
  differ in two documented ways. transitiontrees compares $`G^2`$ with
  the chi-square critical value on $`k - 1`$ degrees of freedom; PST’s
  worked example uses a cutoff equivalent to one degree of freedom.
  transitiontrees computes $`G^2`$ from the empirical child
  distribution, PST from the smoothed one, a difference of the order of
  `ymin` that can change a decision exactly at the cutoff. The suite
  measures the consequence: with PST’s one-degree-of-freedom cutoff, the
  contexts transitiontrees keeps are always a subset of those PST keeps,
  and the two sets are identical when there are two states; with the
  cutoff matched to $`k - 1`$ degrees of freedom, the kept sets are
  identical in at least 85% of configurations, the remainder being
  decisions at the cutoff.
- **Distances, tuning and generation.**
  [`tree_distance()`](https://pak.dynasite.org/transitiontrees/reference/tree_distance.md)
  is a count-weighted symmetric Kullback-Leibler divergence, a different
  quantity from PST’s distance.
  [`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
  uses sequence-level k-fold cross-validation with its own fold
  assignment.
  [`generate_sequences()`](https://pak.dynasite.org/transitiontrees/reference/generate_sequences.md)
  and [`simulate()`](https://rdrr.io/r/stats/simulate.html) consume
  random numbers differently from PST, so they can be compared in
  distribution only.
- **Smoothers.** The `laplace`, `kneser_ney`, `witten_bell` and
  `jelinek_mercer` smoothers have no PST counterpart. The `floor`
  smoother matches PST by default;
  `smoothing = list("floor", ymin = 0.001, rule = "cap")` selects the
  package’s earlier clamp-and-renormalise rule, which differs from PST
  only at contexts with a zero count, by an amount of the order of
  `ymin`.
- **Reshaping.**
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  produces the same set of sequences as
  [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html);
  their order can differ because the two functions sort session
  identifiers differently. Row order does not affect the fit.

## Running the suite

From the root of the source repository, with PST and the other reference
packages available, the suite runs with

``` sh
TRANSITREES_EQUIV_TESTS=true Rscript local_testing_and_equivalence/run-equiv.R
```

which installs the reference packages that are missing (PST from the
CRAN archive), installs the package from source and prints one row per
test with its number of expectations, failures and errors.

## References

Bühlmann, P., & Wyner, A. J. (1999). Variable length Markov chains. *The
Annals of Statistics, 27*(2), 480-513.
<https://doi.org/10.1214/aos/1018031204>

Gabadinho, A., & Ritschard, G. (2016). Analyzing state sequences with
probabilistic suffix trees: The PST R package. *Journal of Statistical
Software, 72*(3). <https://doi.org/10.18637/jss.v072.i03>

Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
probabilistic automata with variable memory length. *Machine Learning,
25*(2-3), 117-149. <https://doi.org/10.1023/A:1026490906255>
