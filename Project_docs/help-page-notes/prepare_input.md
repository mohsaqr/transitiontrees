# Implementation notes removed from `?prepare_input` (2026-11-17)

Still true of R/prepare_input.R; kept for maintainers.

- Base-R reimplementation of the standard long -> sequence reshaping
  (dplyr/tidyr replaced by `order()`/`ave()`/matrix indexing). The
  time-based path is byte-identical to the reference reshaper including
  row order; the no-time path gives the same set of sequences but row
  order can differ (string sort of session ids, so `"session10"` precedes
  `"session2"`). See PARITY.md.
- **Time parsing** (`.pt_parse_time()`): POSIXct and Date pass through;
  numeric (or `is_unix_time = TRUE`) is divided by 1, 1e3 or 1e6 and read
  from 1970-01-01 UTC; character uses `format` if given, else
  `as.POSIXct(x, tz = "UTC")`; any unparsed non-missing value is an error.
- **Session ids**: with `time`, `"<actor> session<k>"` (or `"session<k>"`
  without `actor`); a new session starts at a new actor, a missing gap, or
  a gap > `time_threshold` (strictly greater). With `session`,
  `"<actor> | <session>"` and ordering is `order(actor, session, order)`;
  time is ignored. Without either, the actor id.
- Multiple `actor` columns are pasted with `"-"`.
- `meta` columns: first event's value per session (after sorting), returned
  as a data.frame on `attr(out, "meta")`, aligned to the rows; read by
  `context_tree()` for long-mode `group`/`block`.
- A non-data.frame or zero-row `data` fails a named `stopifnot()`
  (unclassed).
- `.ct_maybe_reshape()` is shared by `context_tree()` and `tune_tree()`;
  naming `actor`/`time`/`order`/`session` without `action` raises
  `transitiontrees_bad_input`.
