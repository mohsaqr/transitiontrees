# Reshape Long Event Data into a Wide Sequence Frame

Turns a long, one-row-per-event table into the wide, one-row-per-
sequence character frame that
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
consumes. Events are grouped by `actor` and ordered by `time` (or
`order`); when `time` is given, each actor's events are split into
**sessions** whenever the gap to the previous event exceeds
`time_threshold` seconds. This mirrors the standard timestamp /
session-splitting rule, in pure base R.

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

  A long-format `data.frame`, one row per event.

- actor:

  Character. Column(s) naming the unit each sequence belongs to (e.g. a
  user id). Several columns are combined with `"-"`. `NULL` (default)
  treats the whole table as one actor.

- time:

  Character. Column holding the event timestamp, used both to order
  events and to split sessions by `time_threshold`. Numeric values are
  read as Unix time. `NULL` (default) uses `order` (or row order)
  instead and does no session splitting.

- action:

  Character. Column holding the event's state / code — the symbol that
  becomes a cell of the sequence. Required.

- order:

  Character. Optional column giving an explicit within- actor ordering
  (used when `time` is absent, or to break time ties). Defaults to row
  order.

- session:

  Character. Optional column giving an explicit session id within an
  actor. If supplied, sessions are taken from it directly and no
  time-gap splitting is done.

- time_threshold:

  Numeric. Seconds; a gap larger than this starts a new session. Default
  `900` (15 minutes), the standard session-splitting threshold.

- format:

  Character. Optional [`strptime`](https://rdrr.io/r/base/strptime.html)
  format for a character `time` column.

- is_unix_time:

  Logical. Force `time` to be read as Unix time. Default `FALSE`.

- unix_time_unit:

  One of `"seconds"` (default), `"milliseconds"`, `"microseconds"`.

- meta:

  Optional character vector of column names to carry through the reshape
  as per-sequence metadata (one value per session, the first event's).
  Returned as a `data.frame` on the `"meta"` attribute of the result,
  aligned to its rows. Used by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  to align column-name `group`/ `block` to the wide rows. Default
  `NULL`.

## Value

A wide character `data.frame`, one row per sequence (session), columns
`T1, T2, ...` holding the ordered states and trailing `NA`s past the end
of each sequence. Row names are the session ids. When `meta` is given,
an aligned per-sequence metadata `data.frame` is attached as
`attr(., "meta")`. Pass it straight to
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## See also

[`context_tree`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)

## Examples

``` r
long <- data.frame(
  user  = c("a","a","a","a","b","b"),
  t     = as.POSIXct("2020-01-01 09:00:00", tz = "UTC") +
            c(0, 60, 3600, 3660, 0, 30),
  state = c("X","Y","X","Z","Y","X"),
  stringsAsFactors = FALSE)
## one-hour gap splits user a into two sessions
wide <- prepare_input(long, actor = "user", time = "t", action = "state")
wide
#>            T1 T2
#> a session1  X  Y
#> a session2  X  Z
#> b session1  Y  X
context_tree(wide, max_depth = 2L, min_count = 1L)
#> <transitiontrees>  3 nodes, depth <= 1, 3 states  [unpruned]
#>   alphabet : X, Y, Z
#>   fit on   : 3 sequences, 6 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 1
#> (start)   n=6      -> X (0.50)
#> |-- X         n=2      -> Y (0.50)
#> `-- Y         n=1      -> X (1.00) 
```
