# Dataset provenance checks for `?trajectories`, `?group_regulation_long`, `?ai_long` and `?engagement` (2026-11-17)

Checks run on the bundled `data/*.rda` files.

- **group_regulation_long**: `identical(group_regulation_long,
  as.data.frame(tna::group_regulation_long))` is `TRUE` (tna 1.3.1; tna
  ships it as a tibble). tna documents it as "generated synthetically from
  `group_regulation`". Also equal to `Nestimate/data/group_regulation_long.rda`
  in the local Nestimate repository. Facts in the help page: 27,533 rows;
  2,000 actors (1,000 High, 1,000 Low; Achiever constant within actor);
  200 groups, one per actor; courses A/B/C (12,390 / 9,626 / 5,517 rows);
  9 actions; 2-26 actions per actor; all timestamps on 2025-01-01. With
  the default `time_threshold = 900`, `prepare_input()` gives 2,000
  sequences (one per actor), max length 26.
- **engagement**: equals `tna::engagement` (an `stslist` with levels
  Active/Average/Disengaged plus `*` and `%`) after converting each column
  to character and setting `*`/`%` to `NA`: `identical()` on the
  dimnames-stripped matrices is `TRUE`, NA positions included. tna's
  source: synthetic, based on Saqr, López-Pernas, Helske & Hrastinski
  (2023), Computers & Education 199, 104787 (DOI verified on Crossref).
  1,000 rows x 25 columns; rows have 23 (145), 24 (155) or 25 (700)
  observed states; NAs only trailing.
- **trajectories**: `identical()` to `Nestimate/data/trajectories.rda`
  (local Nestimate repository). Nestimate's help describes the columns as
  "the 15 successive courses of a programme" but gives no source; no
  source was established, so the help page keeps "Bundled example
  dataset". 138 x 15 character matrix, columns "1".."15"; NAs only
  trailing; 2 all-NA rows; 118 complete rows; 1,870 observed states.
- **ai_long**: `identical()` to `Nestimate/data/ai_long.rda` (local
  Nestimate repository). Nestimate documents it as the AI turns of a
  human-AI vibe-coding study (source given there as a 2026 blog post by
  M. Saqr; not verified here, so not cited). Facts in the help page,
  computed from the data: 8,551 rows; 34 projects; 428 sessions; 5,970
  messages; dates 2026-02-05 to 2026-03-08; timestamps all at 22:00:00
  UTC (midnight UTC+2); 8 codes, 3 clusters; 1-138 rows per session.
