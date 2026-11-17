# Data input: wide sequences, event logs and network objects

A prediction suffix tree is estimated from ordered categorical
sequences: for every position of every sequence, the state that occurs
and the states that precede it. The same sequences reach an analysis in
several shapes. A panel study stores one row per person and one column
per time point; a learning platform or a coding scheme produces an event
log with one row per event, a person, a time stamp and a code; and a
sequence or network analysis that is already under way holds the
sequences inside a fitted object.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
accepts each of these shapes and reduces them to the same internal list
of sequences, so that the tree does not depend on the route by which the
data arrived.

This guide describes each input route, the rules by which an event log
is cut into sequences, the treatment of missing values, and the hand-off
from objects of the `tna`, `Nestimate` and `TraMineR` packages. Each
section confirms that two routes to the same sequences give the same
tree. [A complete
workflow](https://pak.dynasite.org/transitiontrees/articles/workflow.md)
shows the analysis that follows.

``` r

library(transitiontrees)
```

## Data

The package bundles four example datasets, two in each shape.

- **`trajectories`** — a wide character matrix: 138 learners (rows) by
  15 time steps (columns), with the states `Active`, `Average` and
  `Disengaged`; trailing `NA`s mark the end of a learner’s record.
- **`engagement`** — a wide data frame: 1,000 students by 25 weeks (`T1`
  to `T25`) over the same three states; synthetic data from the `tna`
  package.
- **`group_regulation_long`** — an event log of 27,533 collaborative
  regulation moves: the student (`Actor`), the time stamp (`Time`,
  POSIXct), the move (`Action`, nine codes), and the covariates
  `Achiever`, `Group` and `Course`; synthetic data from the `tna`
  package.
- **`ai_long`** — an event log of 8,551 coded moves from projects in
  which people worked with an AI assistant: the project (`project`), an
  explicit session identifier (`session_id`), a Unix time stamp
  (`timestamp`), the move (`code`, eight codes), the order of the move
  within its session (`order_in_session`), and a coarser grouping of the
  codes (`cluster`).

## Wide sequences

In the wide shape each row is one sequence and each column one position,
read from left to right. A character matrix and a data frame are read in
the same way.

``` r

head(engagement, 3)
#>        T1     T2     T3     T4     T5     T6      T7     T8     T9    T10    T11
#> 1 Average Active Active Active Active Active  Active Active Active Active Active
#> 2 Average Active Active Active Active Active  Active Active Active Active Active
#> 3  Active Active Active Active Active Active Average Active Active Active Active
#>          T12        T13        T14        T15        T16        T17        T18        T19
#> 1     Active    Average    Average     Active     Active     Active Disengaged Disengaged
#> 2     Active    Average    Average     Active    Average     Active     Active     Active
#> 3 Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#>          T20     T21        T22        T23     T24  T25
#> 1 Disengaged  Active     Active     Active    <NA> <NA>
#> 2     Active  Active Disengaged Disengaged    <NA> <NA>
#> 3    Average Average Disengaged Disengaged Average <NA>
from_frame <- context_tree(engagement, max_depth = 2, min_count = 5)
from_matrix <- context_tree(as.matrix(engagement), max_depth = 2, min_count = 5)
from_frame
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> `-- Disengaged  n=6598   -> Disengaged (0.79)
#>     |-- Active    n=581    -> Disengaged (0.82)
#>     |-- Average   n=719    -> Disengaged (0.73)
#>     `-- Disengaged  n=4960   -> Disengaged (0.79)
identical(as.data.frame(from_frame), as.data.frame(from_matrix))
#> [1] TRUE
```

The tree is fitted on 1,000 sequences and 24,555 observed states, and
the two containers give identical pathway tables. A sequence may end
early: the `NA` cells after the last observed state are dropped, so
sequences of different lengths need no padding rule. In `trajectories`,
two of the 138 rows contain no observed state at all and do not form a
sequence.

``` r

context_tree(trajectories, max_depth = 1, min_count = 5)
#> <transitiontrees>  4 nodes, depth <= 1, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |-- Average   n=751    -> Average (0.61)
#> `-- Disengaged  n=325    -> Disengaged (0.48)
```

The banner reports 136 sequences and 1,870 observations.

Numeric matrices are refused. In the sibling network packages a square
numeric matrix denotes a transition matrix, and a sequence coded by
integers must be converted to character states explicitly, so that the
alphabet is never inferred from numbers by accident. The refusal is a
classed error, `transitiontrees_bad_input`.

``` r

context_tree(matrix(1:6, nrow = 2))
#> Error:
#> ! 'data' must be a wide data.frame, character/logical matrix, list of character vectors, or a sequence-bearing transition/network object.
```

## Lists of sequences

A list of character vectors holds one sequence per element. The elements
may differ in length.

``` r

sessions <- list(c("plan", "discuss", "consensus", "plan"),
                 c("monitor", "discuss", "plan"),
                 c("plan", "plan", "consensus", "discuss", "plan"))
context_tree(sessions, max_depth = 1, min_count = 1)
#> <transitiontrees>  5 nodes, depth <= 1, 4 states  [unpruned]
#>   alphabet : consensus, discuss, monitor, plan
#>   fit on   : 3 sequences, 12 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 1
#> (start)   n=12     -> plan (0.50)
#> |-- consensus  n=2      -> discuss (0.50)
#> |-- discuss   n=3      -> plan (0.66)
#> |-- monitor   n=1      -> discuss (1.00)
#> `-- plan      n=3      -> consensus (0.33)
```

The three sequences contribute 12 observed states over four codes.

## Missing values and void codes

In wide and list input four codes mark the absence of an observation:
`NA`, the empty string `""`, and the `TraMineR` codes `"%"` (void, the
position after the end of a sequence) and `"*"` (missing). A cell with
any of them is dropped, the observed states on either side are joined,
and the code never becomes a state of the alphabet. The two frames below
hold the same observed states, one with gap codes and one without.

``` r

with_gaps <- data.frame(T1 = c("A", "B", "A"),
                        T2 = c("%", "B", "*"),
                        T3 = c("B", "", "A"),
                        T4 = c(NA, "%", "%"))
without_gaps <- data.frame(T1 = c("A", "B", "A"),
                           T2 = c("B", "B", "A"))
identical(as.data.frame(context_tree(with_gaps, max_depth = 1, min_count = 1)),
          as.data.frame(context_tree(without_gaps, max_depth = 1, min_count = 1)))
#> [1] TRUE
```

Dropping an internal gap joins the states on either side of it, which
creates a transition that was not observed. When gaps inside sequences
are frequent,
[`impute_sequences()`](https://pak.dynasite.org/transitiontrees/reference/impute_sequences.md)
fills them from a fitted tree instead (see [prediction and
scoring](https://pak.dynasite.org/transitiontrees/articles/prediction-and-scoring.md)).

## Event logs

An event log has one row per event.
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
turns it into the wide shape: the events of each `actor` are ordered by
`time` (or by an `order` column) and are cut into sessions wherever the
gap between consecutive events exceeds `time_threshold` seconds (default
900, that is, 15 minutes). Each session becomes one row; the columns are
`T1`, `T2`, and so on.

``` r

wide <- prepare_input(group_regulation_long,
                      actor = "Actor", time = "Time", action = "Action")
dim(wide)
#> [1] 2000   26
subset(head(wide, 3), select = T1:T8)
#>                    T1        T2      T3        T4    T5        T6   T7        T8
#> 1 session1   cohesion consensus discuss synthesis adapt consensus plan consensus
#> 10 session1   emotion  cohesion discuss synthesis  <NA>      <NA> <NA>      <NA>
#> 100 session1     plan consensus    plan      <NA>  <NA>      <NA> <NA>      <NA>
```

The 27,533 events form 2,000 sessions, one per student, the longest with
26 moves; the first eight positions of three sessions are shown. The row
names identify the student and the session. Naming `action` in
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
performs the same reshaping before fitting, so the two routes give the
same tree.

``` r

from_log <- context_tree(group_regulation_long,
                         actor = "Actor", time = "Time", action = "Action",
                         max_depth = 2, min_count = 10)
from_wide <- context_tree(wide, max_depth = 2, min_count = 10)
identical(as.data.frame(from_log), as.data.frame(from_wide))
#> [1] TRUE
```

### The session threshold

The threshold decides what counts as one sequence, and a transition is
counted only within a session. It should be chosen from the meaning of a
pause in the data rather than from its effect on the results. The number
of sessions under three thresholds:

``` r

vapply(c(300, 900, 3600), \(seconds)
       nrow(prepare_input(group_regulation_long, actor = "Actor",
                          time = "Time", action = "Action",
                          time_threshold = seconds)),
       integer(1))
#> [1] 12305  2000  2000
```

A threshold of five minutes cuts the moves of the 2,000 students into
12,305 sessions; 15 minutes and one hour leave one session per student.

### Covariates of a sequence

A covariate recorded on the events, such as `Achiever`, describes the
sequence as a whole. `prepare_input(meta = )` carries the named columns
through the reshaping, one value per session taken from its first event.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
uses the same mechanism when `group` names a column of the log, and fits
one tree per group over a shared alphabet.

``` r

by_achievement <- context_tree(group_regulation_long,
                               actor = "Actor", time = "Time",
                               action = "Action", group = "Achiever",
                               max_depth = 2, min_count = 10)
summary(by_achievement)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC perplexity
#>   High  1000 13721      74         2  FALSE -22787 592 13721 46758 51214       5.26
#>    Low  1000 13812      77         2  FALSE -22984 616 13812 47199 51840       5.28
```

Each achievement group holds 1,000 sessions. [Comparing
groups](https://pak.dynasite.org/transitiontrees/articles/comparing-groups.md)
continues from such a grouped fit.

### Unix time and explicit sessions

A numeric time column is read as Unix time in seconds; `unix_time_unit`
selects milliseconds or microseconds, and `format` parses a character
time stamp. In `ai_long` the time stamp records the day of a session
only: every event of a day carries the same value. Splitting by time
therefore makes one sequence per project and day, with the events of a
day kept in the order of `order_in_session`.

``` r

by_day <- prepare_input(ai_long, actor = "project", time = "timestamp",
                        order = "order_in_session", action = "code")
dim(by_day)
#> [1] 160 227
```

The log also records the session in which each move was made. Naming it
in `session` takes the sessions directly from that column, without
time-based splitting.

``` r

by_session <- prepare_input(ai_long, actor = "project",
                            session = "session_id",
                            order = "order_in_session", action = "code")
dim(by_session)
#> [1] 428 138
```

The 8,551 moves form 160 project-days but 428 sessions. When the data
record sessions, the recorded session is the better unit: a project-day
joins the last move of one session to the first move of the next. The
tree is fitted from the log in the same way.

``` r

ai_tree <- context_tree(ai_long, actor = "project", session = "session_id",
                        order = "order_in_session", action = "code",
                        max_depth = 2, min_count = 5)
state_frequencies(ai_tree)
#>         state count proportion
#> 1     Execute  3258     0.3810
#> 2 Investigate  2317     0.2710
#> 3        Plan  1620     0.1895
#> 4     Explain   524     0.0613
#> 5    Delegate   295     0.0345
#> 6      Repair   257     0.0301
#> 7      Report   181     0.0212
#> 8         Ask    99     0.0116
```

`Execute` and `Investigate` account for 38% and 27% of the moves; `Ask`
and `Report` for about 1% and 2%.

## Objects from tna

A `tna` model keeps the sequences from which it was estimated.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads them from the model and decodes them through the model’s labels,
so that the tree uses the same states as the network.

``` r

engagement_tna <- tna::tna(engagement)
from_tna <- context_tree(engagement_tna, max_depth = 2, min_count = 5)
identical(as.data.frame(from_frame), as.data.frame(from_tna))
#> [1] TRUE
```

An event log prepared by
[`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html)
can be passed to
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
directly; the wide sequences it holds are used. `prepare_data()` reports
its session rule in messages, which are not shown here. The tree is the
same as the one fitted from the log directly: the distance between the
two trees is zero and their fit statistics agree.

``` r

prepared <- tna::prepare_data(group_regulation_long, actor = "Actor",
                              time = "Time", action = "Action")
from_prepared <- context_tree(prepared, max_depth = 2, min_count = 10)
tree_distance(from_log, from_prepared)
#> [1] 0
all.equal(model_fit(from_log), model_fit(from_prepared))
#> [1] TRUE
```

A grouped `tna` model fits one tree per group.

``` r

from_group_model <- context_tree(tna::group_model(prepared, group = "Achiever"),
                                 max_depth = 2, min_count = 10)
summary(from_group_model)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC perplexity
#>   High  1000 13721      74         2  FALSE -22787 592 13721 46758 51214       5.26
#>    Low  1000 13812      77         2  FALSE -22984 616 13812 47199 51840       5.28
```

The two trees are those of the `group = "Achiever"` fit above.

A network estimated from a transition matrix alone carries no sequences.
The sequences cannot be recovered from aggregated transition weights,
and
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
refuses such an object with a `transitiontrees_bad_input` error rather
than fabricating sequences.

``` r

states <- c("Active", "Average", "Disengaged")
weights <- matrix(c(0.6, 0.3, 0.1,
                    0.2, 0.5, 0.3,
                    0.1, 0.3, 0.6), nrow = 3, byrow = TRUE,
                  dimnames = list(states, states))
context_tree(tna::tna(weights))
#> Error:
#> ! This tna object carries no sequence data anywhere (no usable $data / $sequences / $seqdata / embedded netobject). It looks like a pure graph - an aggregated transition network - and transitiontrees fits on raw sequences, not on aggregated transitions: the original sequences cannot be recovered from edge weights (the same reason numeric transition matrices are rejected). Pass a sequence-bearing object instead, e.g. a wide-format transition/network object that still carries its raw sequences, or the original wide sequence data.frame.
```

## Objects from Nestimate

[`Nestimate::build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.html)
and the other `build_*()` constructors return a `netobject` that also
keeps its sequences and a table of node labels, and a grouped
construction returns a `netobject_group`.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads both in the same way as the `tna` objects. The chunk below runs
when Nestimate is installed; it prints `TRUE` when the two routes give
the same pathway table.

``` r

engagement_net <- Nestimate::build_tna(engagement)
from_nestimate <- context_tree(engagement_net, max_depth = 2, min_count = 5)
identical(as.data.frame(from_frame), as.data.frame(from_nestimate))
#> [1] TRUE
```

## State sequence objects from TraMineR

A `TraMineR` state sequence object (class `stslist`) is a data frame of
states in which the end of a shorter sequence is filled with the void
code `"%"`.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads the object, or the character matrix from
[`as.matrix()`](https://rdrr.io/r/base/matrix.html), and drops the void
and missing codes.

`seqdef()` describes the states it found in messages, which are not
shown.

``` r

engagement_seq <- TraMineR::seqdef(engagement)
from_stslist <- context_tree(engagement_seq, max_depth = 2, min_count = 5)
from_seq_matrix <- context_tree(as.matrix(engagement_seq),
                                max_depth = 2, min_count = 5)
identical(as.data.frame(from_frame), as.data.frame(from_stslist))
#> [1] TRUE
identical(as.data.frame(from_frame), as.data.frame(from_seq_matrix))
#> [1] TRUE
```

Both routes give the tree fitted from the wide data frame, with no
alphabet to declare and no void code to remove by hand.
