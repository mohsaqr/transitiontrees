# ---- Imputation of missing states ----

#' Fill in missing states of incomplete sequences
#'
#' @description
#' \code{impute_sequences()} fills in the missing states of incomplete sequences
#' from a fitted tree. Sequences recorded over time often have gaps: positions
#' at which the state was not observed although the sequence continues
#' afterwards. The tree gives, for each gap, a distribution of the missing state
#' conditional on the states observed before it, which is used to fill the gap
#' with the most likely state or with a random draw.
#'
#' Gaps are filled from left to right. Each gap is predicted from the longest
#' context of the tree that matches the states before it, and a filled state
#' becomes part of the history for the gaps that follow. Only gaps followed by
#' an observed state are filled: missing positions at the end of a sequence mark
#' where it stops, and a sequence with no observed state is returned unchanged.
#' The prediction uses past states only; the states observed after a gap do not
#' inform it.
#'
#' With \code{method = "modal"} the result is deterministic and the filled
#' states are the single most likely ones, which understates the variability of
#' the process; with \code{method = "prob"} the states are drawn from the
#' predicted distribution, and repeated imputations with different seeds reflect
#' the uncertainty of the filled values. It returns \code{newdata} in the same
#' form with its internal gaps filled.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Sequences with gaps, marked by \code{NA} or \code{""}: a
#'   list of character vectors, a character matrix or \code{data.frame} with
#'   one sequence per row, or a single character vector.
#' @param method \code{"modal"} (default) fills each gap with the most
#'   likely state; \code{"prob"} draws it from the predicted distribution.
#' @param seed Random seed for \code{method = "prob"}, set locally without
#'   changing the session's random number stream. \code{NULL} (default)
#'   uses the current stream.
#'
#' @return \code{newdata} in the same form (list, matrix, \code{data.frame}
#'   or vector) with its internal gaps filled.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{newdata} is not a list, matrix, \code{data.frame} or character
#' vector.
#'
#' @seealso \code{\link{predict.transitiontrees}()},
#'   \code{\link{generate_sequences}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' gappy <- list(c("Active", NA, "Active", "Average"),
#'               c("Disengaged", "Disengaged", NA, NA, "Average", NA))
#' impute_sequences(tree, gappy)
#' impute_sequences(tree, gappy, method = "prob", seed = 1)
#' @export
impute_sequences <- function(tree, newdata,
                             method = c("modal", "prob"), seed = NULL) {
  .tt_check_class(tree, "transitiontrees", "tree")
  method <- match.arg(method)
  if (!is.null(seed)) {
    ## seed locally without clobbering the caller's RNG stream
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()))
    }
    set.seed(seed)
  }
  alpha <- tree$alphabet

  fill_one <- function(seq_chr) {
    seq_chr <- as.character(seq_chr)
    present <- !is.na(seq_chr) & nzchar(seq_chr)
    if (all(present)) return(seq_chr)
    last_obs <- if (any(present)) max(which(present)) else 0L
    if (last_obs == 0L) return(seq_chr)   # nothing to condition on
    ## Left-to-right with a true sequential dependency: a filled state
    ## becomes context for the gaps that follow it.
    for (t in seq_len(last_obs)) {
      if (present[t]) next
      hist <- seq_chr[seq_len(t - 1L)]
      hist <- hist[!is.na(hist) & nzchar(hist)]
      p    <- tree$nodes[[.ct_match_context(tree, hist)]]$prob
      seq_chr[t] <- if (method == "modal")
        alpha[which.max(p)] else sample(alpha, 1L, prob = p)
      present[t] <- TRUE
    }
    seq_chr
  }

  if (is.list(newdata) && !is.data.frame(newdata)) {
    return(lapply(newdata, fill_one))
  }
  if (is.matrix(newdata) || is.data.frame(newdata)) {
    m <- as.matrix(newdata); storage.mode(m) <- "character"
    filled <- t(vapply(seq_len(nrow(m)), function(i) fill_one(m[i, ]),
                       character(ncol(m))))
    dimnames(filled) <- dimnames(m)
    if (is.data.frame(newdata))
      return(as.data.frame(filled, stringsAsFactors = FALSE))
    return(filled)
  }
  if (is.character(newdata) && is.null(dim(newdata))) {
    return(fill_one(newdata))
  }
  .tt_abort("'newdata' must be a list, matrix, data.frame, or character vector.")
}
