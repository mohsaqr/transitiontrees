# ---- Prediction and generation ----

#' @noRd
.ct_match_context <- function(tree, history) {
  ## Walks history backwards, returning the deepest surviving context
  ## that matches the tail of `history`. Always falls back to root ("").
  history <- as.character(history)
  L <- length(history)
  if (L == 0L) return(.ROOT)
  ## Loop with early return: the deepest matching context wins, so the
  ## search stops at the first hit instead of building every suffix.
  for (k in seq.int(min(tree$max_depth, L), 1L)) {
    ctx <- paste(history[(L - k + 1L):L], collapse = " -> ")
    if (ctx %in% names(tree$nodes)) return(ctx)
  }
  .ROOT
}

#' Predict the next state after each history
#'
#' @description
#' \code{predict()} predicts the state that follows each history from the
#' longest context of a fitted tree that matches the end of that history.
#' Histories longer than the tree's depth are read only as far back as the tree
#' remembers, and histories whose recent states never formed a context fall back
#' to a shorter suffix or to the root. The prediction is therefore always
#' defined, and it uses as much of the past as the data supported when the tree
#' was fitted.
#'
#' It returns, for each history, the probability of every next state, or the
#' most likely next state with \code{type = "class"}. With \code{tidy = TRUE} it
#' returns one row per history with the history, the context that was matched,
#' the most likely next state and its probability, and the probability of every
#' state.
#'
#' To score observed sequences rather than predict from given histories,
#' use \code{\link{score_positions}()}; to look up a single pathway, use
#' \code{\link{query_pathway}()}; to fill in missing states, use
#' \code{\link{impute_sequences}()}.
#'
#' @param object A \code{transitiontrees} object.
#' @param newdata Histories, oldest state first: a list of character
#'   vectors, a wide \code{data.frame} or matrix with one history per row
#'   (missing or empty cells are dropped), or a single character vector.
#'   A history given as one arrow-notation string (\code{"A -> B"}) is
#'   split into its states.
#' @param type \code{"prob"} (default) for next-state probabilities or
#'   \code{"class"} for the most likely next state.
#' @param tidy Return a \code{data.frame} with one row per history.
#'   Default \code{FALSE}.
#' @param ... Not used.
#'
#' @return With \code{tidy = TRUE}, a \code{data.frame} with one row per
#'   history and columns \code{id} (position in \code{newdata}),
#'   \code{history} (arrow notation; \code{"(start)"} when empty),
#'   \code{context} (the matched context), \code{predicted} (most likely
#'   next state), \code{probability} (its probability) and one column per
#'   state with the next-state probabilities. With \code{tidy = FALSE} and
#'   \code{type = "prob"}, a numeric matrix with one row per history and one
#'   column per state, or a named vector when \code{newdata} is a single
#'   character vector; with \code{type = "class"}, a character vector of
#'   most likely next states.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{newdata} is not a list,
#' data frame, matrix or character vector. A warning of class
#' \code{transitiontrees_unknown_state} when a history contains a state
#' outside the tree's alphabet; the prediction then uses the longest
#' context after the last unknown state.
#'
#' @seealso \code{\link{query_pathway}()}, \code{\link{score_positions}()},
#'   \code{\link{impute_sequences}()}, \code{\link{generate_sequences}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' histories <- list(c("Active", "Active"),
#'                   c("Average", "Disengaged", "Disengaged"),
#'                   character(0))
#' predict(tree, newdata = histories, tidy = TRUE)
#' predict(tree, newdata = histories, type = "class")
#' @export
predict.transitiontrees <- function(object, newdata, type = c("prob", "class"),
                            tidy = FALSE, ...) {
  type <- match.arg(type)
  stopifnot("`tidy` must be TRUE or FALSE" =
              is.logical(tidy) && length(tidy) == 1L && !is.na(tidy))
  ## A bare character vector is the interactive shortcut; a character
  ## *matrix* is a container (it has a dim) and must yield a matrix, so
  ## exclude it here even though is.character() is TRUE for it.
  bare_vec <- is.character(newdata) && !is.list(newdata) &&
    is.null(dim(newdata))

  histories <- if (is.list(newdata) && !is.data.frame(newdata)) {
    lapply(newdata, as.character)
  } else if (is.data.frame(newdata) || is.matrix(newdata)) {
    m <- as.matrix(newdata); storage.mode(m) <- "character"
    lapply(seq_len(nrow(m)), function(i) {
      r <- m[i, ]
      r[!is.na(r) & nzchar(r)]
    })
  } else if (bare_vec) {
    list(newdata)
  } else {
    .tt_abort("'newdata' must be a list, data.frame, matrix, or character vector.")
  }
  ## A history given as one arrow-notation string ("A -> B") is the
  ## package's own pathway notation: split it into states.
  histories <- lapply(histories, function(h)
    if (length(h) == 1L && grepl(" -> ", h, fixed = TRUE))
      strsplit(h, " -> ", fixed = TRUE)[[1L]] else h)
  ## A state outside the alphabet ends the usable history (no context
  ## can contain it); say so instead of silently predicting from a
  ## shorter context.
  unknown <- setdiff(unique(unlist(histories)), object$alphabet)
  if (length(unknown))
    .tt_warn("History state(s) not in the tree's alphabet: ",
             paste(utils::head(unknown, 10L), collapse = ", "),
             ". Prediction uses the longest context after the last unknown ",
             "state.", class = "transitiontrees_unknown_state")

  probs <- t(vapply(histories, function(h) {
    ctx  <- .ct_match_context(object, h)
    info <- object$nodes[[ctx]]
    if (is.null(info)) {
      ## Can only happen if root is missing — return uniform
      rep(1 / length(object$alphabet), length(object$alphabet))
    } else info$prob
  }, numeric(length(object$alphabet))))
  colnames(probs) <- object$alphabet

  if (isTRUE(tidy)) {
    contexts <- vapply(histories, function(h) .ct_match_context(object, h),
                       character(1))
    hist_lab <- vapply(histories, function(h)
      if (length(h)) paste(h, collapse = " -> ") else .ROOT_LABEL,
      character(1))
    top <- max.col(probs, ties.method = "first")
    fixed <- data.frame(
      id          = seq_along(histories),
      history     = hist_lab,
      context     = ifelse(contexts == .ROOT, .ROOT_LABEL, contexts),
      predicted   = object$alphabet[top],
      probability = probs[cbind(seq_len(nrow(probs)), top)],
      stringsAsFactors = FALSE, row.names = NULL)
    out <- cbind(fixed, as.data.frame(probs, stringsAsFactors = FALSE,
                                      optional = TRUE))
    names(out) <- make.unique(names(out))
    rownames(out) <- NULL
    return(out)
  }

  if (type == "class") {
    return(object$alphabet[max.col(probs, ties.method = "first")])
  }
  ## Schema-stable: container inputs always yield a matrix. Only a
  ## bare character vector — an interactive shortcut — collapses to
  ## a named vector.
  if (bare_vec) return(setNames(as.vector(probs), object$alphabet))
  probs
}

#' Generate sequences from a fitted transition tree
#'
#' @description
#' \code{generate_sequences()} draws new sequences from the process described by
#' a fitted tree. Each sequence starts from a state drawn from the overall
#' distribution of states (or from a given first state), and each further state
#' is drawn from the next-state distribution of the longest context of the tree
#' that matches the states generated so far. Generated sequences can be compared
#' with the observed ones to check whether the tree reproduces features of the
#' data that it was not fitted to directly, such as run lengths or the time
#' spent in each state.
#'
#' Draws use the session's random number generator; call
#' \code{\link{set.seed}()} first, or use
#' \code{\link{simulate.transitiontrees}()} with \code{seed}, for reproducible
#' sequences. It returns one sequence per row, as a character matrix or, with
#' \code{tidy = TRUE}, a data frame.
#'
#' @param tree A \code{transitiontrees} object.
#' @param n Number of sequences. Default 5.
#' @param length Number of states in each sequence. Default 10.
#' @param start First state of each sequence, a character vector of length
#'   \code{n}. \code{NULL} (default) draws it from the root distribution.
#' @param tidy Return a \code{data.frame} rather than a matrix. Default
#'   \code{FALSE}.
#'
#' @return With \code{tidy = TRUE}, a \code{data.frame} with one row per
#'   sequence and columns \code{T1} to \code{T<length>}, which
#'   \code{\link{context_tree}()} accepts as input. With
#'   \code{tidy = FALSE}, a character matrix with the same layout.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree,
#' \code{n} or \code{length} is not a single positive number, or
#' \code{start} does not have length \code{n}.
#'
#' @seealso \code{\link{simulate.transitiontrees}()},
#'   \code{\link{predict.transitiontrees}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' set.seed(1)
#' generate_sequences(tree, n = 5, length = 12, tidy = TRUE)
#' @export
generate_sequences <- function(tree, n = 5L, length = 10L,
                               start = NULL, tidy = FALSE) {
  .tt_check_class(tree, "transitiontrees", "tree")
  stopifnot("`tidy` must be TRUE or FALSE" =
              is.logical(tidy) && base::length(tidy) == 1L && !is.na(tidy))
  ## base::length() because the `length` argument shadows the function.
  if (!is.numeric(n) || base::length(n) != 1L || is.na(n) || n < 1)
    .tt_abort("'n' must be a single positive integer.")
  if (!is.numeric(length) || base::length(length) != 1L ||
      is.na(length) || length < 1)
    .tt_abort("'length' must be a single positive integer.")
  n <- as.integer(n); length <- as.integer(length)
  alpha <- tree$alphabet

  draw <- function(p) sample(alpha, 1L, prob = p)

  if (is.null(start)) {
    start <- replicate(n, draw(tree$nodes[[.ROOT]]$prob))
  } else {
    start <- as.character(start)
    if (length(start) != n)
      .tt_abort("'start' must have length 'n'.")
  }

  out <- matrix(NA_character_, n, length)
  out[, 1L] <- start
  ## length == 1 is the start column alone; seq.int(2L, 1L) counts *down*,
  ## so only iterate when there is a second position to fill.
  if (length >= 2L) {
    ## Genuinely sequential: each state is drawn from the context formed by
    ## the states already drawn, and the draw order (position, then
    ## sequence) fixes the RNG stream that seed-based tests rely on.
    for (t in seq.int(2L, length)) {
      for (i in seq_len(n)) {
        hist <- out[i, seq_len(t - 1L)]
        ctx  <- .ct_match_context(tree, hist)
        info <- tree$nodes[[ctx]]
        out[i, t] <- draw(info$prob)
      }
    }
  }
  if (isTRUE(tidy)) {
    colnames(out) <- paste0("T", seq_len(length))
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  out
}

#' Simulate sequences from a fitted transition tree
#'
#' @description
#' \code{simulate()} draws \code{nsim} sequences from a fitted tree as
#' \code{\link{generate_sequences}()} does, and sets the random seed first when
#' \code{seed} is given. It is the \code{\link[stats]{simulate}()} method for
#' fitted trees.
#'
#' @param object A \code{transitiontrees} object.
#' @param nsim Number of sequences. Default 5.
#' @param seed Random seed passed to \code{\link{set.seed}()}, or
#'   \code{NULL} (default) to use the current state of the generator.
#' @param length Number of states in each sequence. Default 10.
#' @param start First state of each sequence; see
#'   \code{\link{generate_sequences}()}. Default \code{NULL}.
#' @param tidy Return a \code{data.frame} rather than a matrix. Default
#'   \code{FALSE}.
#' @param ... Not used.
#'
#' @return As \code{\link{generate_sequences}()}: with \code{tidy = TRUE}, a
#'   \code{data.frame} with one row per sequence and columns \code{T1} to
#'   \code{T<length>}; otherwise a character matrix with the same layout.
#'
#' @section Conditions:
#' As \code{\link{generate_sequences}()}.
#'
#' @seealso \code{\link{generate_sequences}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' simulate(tree, nsim = 3, seed = 42, length = 8, tidy = TRUE)
#' @export
simulate.transitiontrees <- function(object, nsim = 5L, seed = NULL,
                              length = 10L, start = NULL, tidy = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)
  generate_sequences(object, n = nsim, length = length, start = start,
                     tidy = tidy)
}
