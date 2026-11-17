# Reshape a long event log into one sequence per row

`prepare_input()` reshapes a long event log into one sequence per row:
the wide table that
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
fits. Process data are usually recorded as event logs, one row per event
with the actor, a timestamp and the action or state, while sequence
models need the events of each actor, or of each working session, in
temporal order as one sequence.

Events are grouped by `actor` and ordered by `time` (ties broken by
`order`, then by row order). When `time` is given, each actor's events
are split into sessions wherever the gap between consecutive events
exceeds `time_threshold` seconds, so that a sequence covers one
continuous episode of activity rather than all of an actor's events over
days or weeks. A `session` column, when given, defines the sessions
directly, and events within a session are then ordered by `order` or row
order, not by `time`. Without `time` or `session`, each actor
contributes one sequence in the order of `order` or of the rows.

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
calls `prepare_input()` itself when its `action` argument is set, so the
reshaping is needed separately only to inspect the sequences or to pass
them to other functions. The choice of `time_threshold` defines what
counts as a sequence and changes every result fitted on it; it is worth
checking that conclusions hold for neighboring thresholds.

## Usage

``` r
prepare_input(
  data,
  actor = NULL,
  time = NULL,
  action = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  format = NULL,
  is_unix_time = FALSE,
  unix_time_unit = c("seconds", "milliseconds", "microseconds"),
  meta = NULL
)
```

## Arguments

- data:

  A `data.frame` with one row per event.

- actor:

  Name of the column, or columns, identifying whose events form a
  sequence. Several columns are combined. `NULL` (default) treats all
  rows as one actor.

- time:

  Name of the timestamp column: date-times, dates, character timestamps,
  or numbers read as Unix time. `NULL` (default) orders by `order` or
  row order and does not split sessions.

- action:

  Name of the column holding the state of each event. Required.

- order:

  Name of a column giving the order of events within an actor, used
  without `time` or to break ties in time. `NULL` (default) uses row
  order.

- session:

  Name of a column identifying sessions within an actor. When given,
  sessions are taken from it and `time_threshold` is not used. Default
  `NULL`.

- time_threshold:

  Gap in seconds between consecutive events that starts a new session.
  Default 900 (15 minutes).

- format:

  [`strptime`](https://rdrr.io/r/base/strptime.html) format for a
  character `time` column. `NULL` (default) accepts ISO 8601 timestamps.

- is_unix_time:

  Read `time` as Unix time. Default `FALSE`; numeric columns are read as
  Unix time in any case.

- unix_time_unit:

  Unit of Unix time: `"seconds"` (default), `"milliseconds"` or
  `"microseconds"`.

- meta:

  Names of columns to carry over as one value per sequence (the value of
  its first event), used by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  to align `group` and `block` with the sequences. Default `NULL`.

## Value

A character `data.frame` with one row per sequence, named by actor and
session, and columns `T1`, `T2`, ... holding the states in order;
sequences shorter than the longest are padded with `NA`. It can be
passed directly to
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Conditions

`transitiontrees_bad_input` when `action` or another named column is not
in `data`, `time_threshold` is not a single non-negative number, or
`time` cannot be parsed.

## See also

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
which accepts the same arguments.

## Examples

``` r
sessions <- prepare_input(group_regulation_long, actor = "Actor",
                          time = "Time", action = "Action")
dim(sessions)
#> [1] 2000   26
head(sessions, 3)
#>                    T1        T2      T3        T4    T5        T6   T7
#> 1 session1   cohesion consensus discuss synthesis adapt consensus plan
#> 10 session1   emotion  cohesion discuss synthesis  <NA>      <NA> <NA>
#> 100 session1     plan consensus    plan      <NA>  <NA>      <NA> <NA>
#>                     T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19
#> 1 session1   consensus <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> 10 session1       <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> 100 session1      <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#>               T20  T21  T22  T23  T24  T25  T26
#> 1 session1   <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> 10 session1  <NA> <NA> <NA> <NA> <NA> <NA> <NA>
#> 100 session1 <NA> <NA> <NA> <NA> <NA> <NA> <NA>

# Explicit session identifiers and within-session order
turns <- prepare_input(ai_long, actor = "session_id",
                       order = "order_in_session", action = "code")
dim(turns)
#> [1] 428 138
```
