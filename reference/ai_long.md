# AI turns in human-AI programming sessions (long event log)

`ai_long` is a coded dataset of the turns of an AI coding assistant in
428 human-AI pair-programming sessions across 34 projects, held between
5 February and 8 March 2026, one row per code. A turn (message) can
carry several codes, in the order given by `code_order`. Each code is
one of eight assistant actions (`Ask`, `Delegate`, `Execute`, `Explain`,
`Investigate`, `Plan`, `Repair`, `Report`), grouped into three clusters
(`Action`, `Communication`, `Repair`).

Timestamps have day resolution (midnight, UTC+2, of the session date),
so the order of events within a day comes from `order_in_session`, which
counts the human turns of the session as well and therefore has gaps.
Use `actor = "session_id"` with `order = "order_in_session"` for one
sequence per session, or `time = "timestamp"` to split sequences at
changes of day.

## Usage

``` r
ai_long
```

## Format

A `data.frame` with 8,551 rows (codes) and 9 columns:

- message_id:

  Integer identifier of the turn; the codes of one turn share it.

- project:

  Project identifier, `"Project_1"` to `"Project_34"`.

- session_id:

  Session identifier, a 12-character hash.

- timestamp:

  Integer Unix time in seconds (day resolution).

- session_date:

  Date of the row, `"YYYY-MM-DD"`.

- code:

  The assistant action, the state of the sequence.

- cluster:

  The cluster of the code.

- code_order:

  Order of the code within its turn, from 1.

- order_in_session:

  Position of the row in its session, counting human and AI rows
  together.

## Source

Bundled example dataset. The same data are distributed with the
Nestimate package as `ai_long`.

## See also

[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md),
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Examples

``` r
tree <- context_tree(ai_long, actor = "session_id",
                     order = "order_in_session", action = "code",
                     max_depth = 2)
state_frequencies(tree)
#>         state count proportion
#> 1     Execute  3258 0.38100807
#> 2 Investigate  2317 0.27096246
#> 3        Plan  1620 0.18945153
#> 4     Explain   524 0.06127938
#> 5    Delegate   295 0.03449889
#> 6      Repair   257 0.03005496
#> 7      Report   181 0.02116711
#> 8         Ask    99 0.01157759
common_pathways(tree, top = 5, depth = 2)
#>                  pathway depth count likely_next next_probability  divergence
#> 1     Execute -> Execute     2  1503     Execute        0.5522289 0.006393297
#> 2    Investigate -> Plan     2   921     Execute        0.4668838 0.004718142
#> 3 Execute -> Investigate     2   727        Plan        0.3438790 0.031815394
#> 4        Plan -> Execute     2   674     Execute        0.4569733 0.016165977
#> 5 Investigate -> Execute     2   532     Execute        0.5112782 0.007986774
#>   changes_prediction
#> 1              FALSE
#> 2              FALSE
#> 3              FALSE
#> 4              FALSE
#> 5              FALSE
```
