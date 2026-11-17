# ---- prepare_input(): long-format -> wide sequence frame ----
#
# Base-R reimplementation of standard long->sequence reshaping.
# transitiontrees stays pure
# base R, so the usual dplyr/tidyr machinery is replaced by
# order() / ave() / matrix-indexing, but the session-splitting logic is
# identical: within each actor, a new session begins when the gap to the
# previous timestamp exceeds `time_threshold` seconds.

#' Parse a time column to POSIXct.
#' @noRd
.pt_parse_time <- function(x, format = NULL, is_unix_time = FALSE,
                           unix_time_unit = "seconds") {
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.POSIXct(x))
  if (inherits(x, "Date"))                 return(as.POSIXct(x))
  if (is.numeric(x) || isTRUE(is_unix_time)) {
    div <- switch(unix_time_unit,
                  seconds = 1, milliseconds = 1e3, microseconds = 1e6,
                  .tt_abort("'unix_time_unit' must be 'seconds', ",
                       "'milliseconds', or 'microseconds'."))
    return(as.POSIXct(as.numeric(x) / div, origin = "1970-01-01",
                      tz = "UTC"))
  }
  x <- as.character(x)
  if (!is.null(format))
    return(as.POSIXct(x, format = format, tz = "UTC"))
  parsed <- as.POSIXct(x, tz = "UTC")          # ISO-8601 / standard
  if (all(is.na(parsed) == is.na(x)) == FALSE || anyNA(parsed[!is.na(x)]))
    .tt_abort("Could not parse 'time'; pass an explicit 'format' ",
         "(see strptime) or 'is_unix_time = TRUE'.")
  parsed
}

#' Reshape a long event log into one sequence per row
#'
#' @description
#' \code{prepare_input()} reshapes a long event log into one sequence per row:
#' the wide table that \code{\link{context_tree}()} fits. Process data are
#' usually recorded as event logs, one row per event with the actor, a timestamp
#' and the action or state, while sequence models need the events of each actor,
#' or of each working session, in temporal order as one sequence.
#'
#' Events are grouped by \code{actor} and ordered by \code{time} (ties
#' broken by \code{order}, then by row order). When \code{time} is given,
#' each actor's events are split into sessions wherever the gap between
#' consecutive events exceeds \code{time_threshold} seconds, so that a
#' sequence covers one continuous episode of activity rather than all of an
#' actor's events over days or weeks. A \code{session} column, when given,
#' defines the sessions directly, and events within a session are then
#' ordered by \code{order} or row order, not by \code{time}. Without \code{time} or \code{session},
#' each actor contributes one sequence in the order of \code{order} or of
#' the rows.
#'
#' \code{context_tree()} calls \code{prepare_input()} itself when its
#' \code{action} argument is set, so the reshaping is needed separately only
#' to inspect the sequences or to pass them to other functions. The choice
#' of \code{time_threshold} defines what counts as a sequence and changes
#' every result fitted on it; it is worth checking that conclusions hold
#' for neighboring thresholds.
#'
#' @param data A \code{data.frame} with one row per event.
#' @param actor Name of the column, or columns, identifying whose events
#'   form a sequence. Several columns are combined. \code{NULL} (default)
#'   treats all rows as one actor.
#' @param time Name of the timestamp column: date-times, dates, character
#'   timestamps, or numbers read as Unix time. \code{NULL} (default) orders
#'   by \code{order} or row order and does not split sessions.
#' @param action Name of the column holding the state of each event.
#'   Required.
#' @param order Name of a column giving the order of events within an
#'   actor, used without \code{time} or to break ties in time. \code{NULL}
#'   (default) uses row order.
#' @param session Name of a column identifying sessions within an actor.
#'   When given, sessions are taken from it and \code{time_threshold} is not
#'   used. Default \code{NULL}.
#' @param time_threshold Gap in seconds between consecutive events that
#'   starts a new session. Default 900 (15 minutes).
#' @param format \code{\link{strptime}} format for a character \code{time}
#'   column. \code{NULL} (default) accepts ISO 8601 timestamps.
#' @param is_unix_time Read \code{time} as Unix time. Default \code{FALSE};
#'   numeric columns are read as Unix time in any case.
#' @param unix_time_unit Unit of Unix time: \code{"seconds"} (default),
#'   \code{"milliseconds"} or \code{"microseconds"}.
#' @param meta Names of columns to carry over as one value per sequence (the
#'   value of its first event), used by \code{\link{context_tree}()} to
#'   align \code{group} and \code{block} with the sequences. Default
#'   \code{NULL}.
#'
#' @return A character \code{data.frame} with one row per sequence, named by
#'   actor and session, and columns \code{T1}, \code{T2}, \ldots holding the
#'   states in order; sequences shorter than the longest are padded with
#'   \code{NA}. It can be passed directly to \code{\link{context_tree}()}.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{action} or another named
#' column is not in \code{data}, \code{time_threshold} is not a single
#' non-negative number, or \code{time} cannot be parsed.
#'
#' @seealso \code{\link{context_tree}()}, which accepts the same arguments.
#'
#' @examples
#' sessions <- prepare_input(group_regulation_long, actor = "Actor",
#'                           time = "Time", action = "Action")
#' dim(sessions)
#' head(sessions, 3)
#'
#' # Explicit session identifiers and within-session order
#' turns <- prepare_input(ai_long, actor = "session_id",
#'                        order = "order_in_session", action = "code")
#' dim(turns)
#' @export
prepare_input <- function(data, actor = NULL, time = NULL, action = NULL,
                         order = NULL, session = NULL,
                         time_threshold = 900,
                         format = NULL, is_unix_time = FALSE,
                         unix_time_unit = c("seconds", "milliseconds",
                                            "microseconds"),
                         meta = NULL) {
  stopifnot("`data` must be a data.frame with at least one row" =
              is.data.frame(data) && nrow(data) > 0L)
  unix_time_unit <- match.arg(unix_time_unit)
  if (is.null(action) || !action %in% names(data))
    .tt_abort("'action' must name a column of 'data'.")
  missing_cols <- setdiff(c(actor, time, order, session, meta), names(data))
  if (length(missing_cols))
    .tt_abort("'", paste(missing_cols, collapse = "', '"),
         "' is not a column of 'data'.")
  if (!is.numeric(time_threshold) || length(time_threshold) != 1L ||
      time_threshold < 0)
    .tt_abort("'time_threshold' must be a single non-negative number.")

  n   <- nrow(data)
  act <- as.character(data[[action]])

  ## actor identifier (default: a single actor)
  default_actor <- is.null(actor)
  actor_vec <- if (default_actor) rep("session", n)
    else if (length(actor) > 1L)
      do.call(paste, c(lapply(actor, function(a) data[[a]]), sep = "-"))
    else as.character(data[[actor]])

  ord <- if (is.null(order)) seq_len(n) else data[[order]]

  if (!is.null(session)) {
    ## Explicit session column: order within
    ## actor+session, take the session id directly.
    sess <- paste0(actor_vec, " | ", as.character(data[[session]]))
    o    <- order(actor_vec, sess, ord)
    session_id <- sess[o]
    act_s      <- act[o]
  } else if (!is.null(time)) {
    ## Timestamp logic: sort within actor by (time, order), then split
    ## into sessions on gaps > time_threshold (the standard rule).
    tt <- .pt_parse_time(data[[time]], format, is_unix_time, unix_time_unit)
    o  <- order(actor_vec, tt, ord)
    actor_s <- actor_vec[o]; tt_s <- tt[o]; act_s <- act[o]
    gap <- c(NA_real_, as.numeric(difftime(tt_s[-1L], tt_s[-n],
                                           units = "secs")))
    new_actor   <- c(TRUE, actor_s[-1L] != actor_s[-n])
    new_session <- new_actor | is.na(gap) | gap > time_threshold
    session_nr  <- stats::ave(as.integer(new_session), actor_s,
                              FUN = cumsum)
    session_id  <- if (default_actor) paste0("session", session_nr)
      else paste0(actor_s, " session", session_nr)
  } else {
    ## No time and no session: one sequence per actor, ordered by `order`.
    o <- order(actor_vec, ord)
    session_id <- actor_vec[o]
    act_s      <- act[o]
  }

  ## Position within each session, then scatter into a wide matrix.
  seqpos   <- stats::ave(seq_along(session_id), session_id,
                         FUN = seq_along)
  ## rows sorted by session id as a string, matching the standard
  ## reshaping order (so "session10" precedes "session2").
  sessions <- sort(unique(session_id))
  max_len  <- max(seqpos)
  wide <- matrix(NA_character_, length(sessions), max_len,
                 dimnames = list(sessions, paste0("T", seq_len(max_len))))
  wide[cbind(match(session_id, sessions), seqpos)] <- act_s
  out <- as.data.frame(wide, stringsAsFactors = FALSE)

  ## Carry per-sequence metadata through the reshape: one value per
  ## session (first event's), aligned to the wide rows -- so callers get
  ## group/block aligned by construction instead of re-deriving it.
  if (!is.null(meta)) {
    first  <- !duplicated(session_id)
    md     <- data[o, meta, drop = FALSE][first, , drop = FALSE]
    attr(out, "meta") <- md[match(sessions, session_id[first]), ,
                            drop = FALSE]
  }
  out
}

#' Reshape long-format \code{data} to a wide sequence frame when
#' \code{action} is named, else pass it through. Shared by the functions
#' that accept raw data (\code{context_tree}, \code{tune_tree},
#' \code{compare_smoothing}) so they all gain long-format input
#' identically. Naming \code{actor}/\code{time}/\code{order}/\code{session}
#' without \code{action} is an error (the state column must be named).
#' @noRd
.ct_maybe_reshape <- function(data, actor, time, action, order, session,
                              time_threshold, meta = NULL) {
  if (!is.null(action))
    return(prepare_input(data, actor = actor, time = time, action = action,
                         order = order, session = session,
                         time_threshold = time_threshold, meta = meta))
  if (!is.null(actor) || !is.null(time) || !is.null(session) ||
      !is.null(order))
    .tt_abort("To fit from long-format data, name the state column via ",
         "'action ='.")
  data
}
