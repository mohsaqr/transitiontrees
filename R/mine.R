# ---- Mining contexts and sequences ----

#' Find the contexts that make a given next state likely or unlikely
#'
#' @description
#' \code{mine_contexts()} finds the contexts that make a given next state likely
#' or unlikely. It scans every context of a fitted tree for the probability it
#' gives one chosen next state and returns the contexts whose probability lies
#' in a requested range, from the highest probability down. This serves
#' questions about one outcome state, such as after which histories
#' disengagement is most likely, or after which histories a learner almost never
#' returns to planning.
#'
#' Each row reports the context, its depth and count, the probability of
#' the chosen state after it, and whether that state is the context's most
#' likely next state. Probabilities in rare contexts are imprecise; set
#' \code{min_count} to keep only contexts with enough occurrences. To
#' compare the state's probability across groups, fit a grouped tree and
#' use \code{\link{compare_groups}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param state The next state of interest, one of the tree's states.
#' @param min_prob,max_prob Keep contexts whose probability of \code{state}
#'   is at least \code{min_prob} and at most \code{max_prob}. \code{NULL}
#'   (default) leaves that bound open.
#' @param min_count Keep contexts that occurred at least this many times.
#'   Default 1.
#'
#' @return A \code{data.frame} with one row per context and columns
#'   \code{pathway}, \code{depth}, \code{count}, \code{state}, \code{prob}
#'   (probability of \code{state} after the context) and \code{is_modal}
#'   (whether \code{state} is the most likely next state), sorted by
#'   \code{prob} and then \code{count}, both decreasing. No matching context
#'   gives a zero-row table with these columns.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{state} is not a single one of its states.
#'
#' @seealso \code{\link{query_pathway}()}, \code{\link{sharp_pathways}()},
#'   \code{\link{mine_sequences}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' mine_contexts(tree, state = "Disengaged", min_prob = 0.4, min_count = 10)
#' mine_contexts(tree, state = "Active", max_prob = 0.05)
#' @export
mine_contexts <- function(tree, state, min_prob = NULL, max_prob = NULL,
                          min_count = 1L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  alpha <- tree$alphabet
  if (missing(state) || length(state) != 1L || !state %in% alpha)
    .tt_abort("'state' must be a single state from the alphabet: ",
         paste(alpha, collapse = ", "))
  sidx      <- match(state, alpha)
  min_count <- as.integer(min_count)

  empty <- data.frame(pathway = character(0), depth = integer(0),
                      count = numeric(0), state = character(0),
                      prob = numeric(0), is_modal = logical(0),
                      stringsAsFactors = FALSE)

  rows <- lapply(names(tree$nodes), function(ctx) {
    nd <- tree$nodes[[ctx]]
    if (nd$n < min_count) return(NULL)
    data.frame(
      pathway  = if (identical(ctx, .ROOT)) .ROOT_LABEL else ctx,
      depth    = nd$depth,
      count    = as.numeric(nd$n),   # numeric, matching the empty template
      state    = state,
      prob     = nd$prob[sidx],
      is_modal = which.max(nd$prob) == sidx,
      stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(empty)

  out <- do.call(rbind, rows)
  if (!is.null(min_prob)) out <- out[out$prob >= min_prob, , drop = FALSE]
  if (!is.null(max_prob)) out <- out[out$prob <= max_prob, , drop = FALSE]
  out <- out[order(-out$prob, -out$count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Find the least and most typical sequences
#'
#' @description
#' \code{mine_sequences()} ranks sequences by how well a fitted tree predicts
#' them and returns the top \code{n}. With \code{which = "surprising"}, the
#' sequences with the highest perplexity come first: those whose transitions the
#' tree finds least expected, which are candidates for atypical cases or for
#' patterns that the fitting data did not contain. With \code{which =
#' "expected"}, the sequences with the lowest perplexity come first: the most
#' typical sequences under the tree.
#'
#' The ranking is by per-sequence perplexity from
#' \code{\link{score_sequences}()}, which adjusts for sequence length. When
#' \code{newdata} are the fitting sequences themselves, the ranking
#' describes typicality within the data; with other sequences it describes
#' how well the tree generalizes to them. To find where within a sequence
#' the unexpected moves occur, use \code{\link{score_positions}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Sequences to rank, in any form accepted by
#'   \code{\link{context_tree}()}.
#' @param n Number of sequences to return. Default 10.
#' @param which \code{"surprising"} (default; highest perplexity first) or
#'   \code{"expected"} (lowest perplexity first).
#'
#' @return A \code{data.frame} with the columns of
#'   \code{\link{score_sequences}()} (\code{sequence_id}, \code{n_scored},
#'   \code{log_lik}, \code{perplexity}), one row per sequence, at most
#'   \code{n} rows.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{newdata} has no usable sequence.
#'
#' @seealso \code{\link{score_sequences}()}, \code{\link{score_positions}()},
#'   \code{\link{mine_contexts}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' mine_sequences(tree, engagement, n = 5)
#' mine_sequences(tree, engagement, n = 5, which = "expected")
#' @export
mine_sequences <- function(tree, newdata, n = 10L,
                           which = c("surprising", "expected")) {
  .tt_check_class(tree, "transitiontrees", "tree")
  which <- match.arg(which)
  sc <- score_sequences(tree, newdata)
  if (nrow(sc) == 0L) return(sc)
  ord <- if (which == "surprising") order(-sc$perplexity) else
    order(sc$perplexity)
  out <- utils::head(sc[ord, , drop = FALSE], as.integer(n))
  rownames(out) <- NULL
  out
}
