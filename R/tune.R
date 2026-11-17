# ---- tune_tree(): k-fold cross-validated hyperparameter selection ----

#' @noRd
.pt_make_folds <- function(n, k, seed = 1L) {
  set.seed(seed)
  ord  <- sample.int(n)
  fold <- ((seq_len(n) - 1L) %% k) + 1L
  split(ord, fold[seq_len(n)])
}

#' @noRd
.pt_fit_score_fold <- function(train_seqs, test_seqs, alphabet, cfg) {
  tr <- context_tree(train_seqs,
                     max_depth = cfg$max_depth,
                     min_count = cfg$min_count,
                     smoothing = cfg$smoothing,
                     alphabet  = alphabet)
  if (isTRUE(cfg$prune))
    tr <- prune_tree(tr, criterion = "G2", alpha = cfg$alpha)
  ll <- logLik.transitiontrees(tr, newdata = test_seqs)
  list(ll = as.numeric(ll), n = attr(ll, "nobs"),
       n_nodes = length(tr$nodes))
}

#' @noRd
.pt_smoothing_grid <- function(spec) {
  ## Normalise a user-supplied smoothing grid to a list of resolved
  ## smoothing specs. Accepts a character vector of method names, or a
  ## list of method-name-or-list elements.
  if (is.character(spec)) return(lapply(spec, .pt_resolve_smoothing))
  if (is.list(spec)) {
    if (!is.null(spec$method) || (length(spec) >= 1L &&
        is.character(spec[[1L]]) && length(spec[[1L]]) == 1L &&
        spec[[1L]] %in% names(.pt_smoothing_defaults)))
      return(list(.pt_resolve_smoothing(spec)))
    return(lapply(spec, .pt_resolve_smoothing))
  }
  .tt_abort("'smoothing' must be a character vector or a list of specs.")
}

#' @noRd
.pt_smoothing_label <- function(sm) {
  extras <- sm[setdiff(names(sm), "method")]
  if (length(extras) == 0L) return(sm$method)
  paste0(sm$method, "(",
         paste(names(extras), unlist(extras), sep = "=",
               collapse = ", "),
         ")")
}

#' Select context tree settings by cross-validated perplexity
#'
#' @description
#' \code{tune_tree()} selects the settings of a variable-order context tree by
#' cross-validated perplexity. Four settings decide how much history the tree
#' uses and how far it trusts sparse histories: the maximum depth, the minimum
#' count a context needs to be kept, the smoothing of the next-state
#' distributions, and whether the tree is pruned. A deeper, less smoothed tree
#' reproduces the sequences it was fitted to more closely, but it does not
#' necessarily predict new sequences better. Held-out perplexity measures the
#' second property: it is the exponential of the mean negative log-likelihood
#' per scored transition in sequences the model has not seen, and can be read as
#' the effective number of next states the model is choosing between.
#'
#' It fits a tree for every combination of the supplied values and scores it by
#' k-fold cross-validation. Folds are formed from whole sequences, so a held-out
#' sequence contributes no transitions to the fit, and every configuration is
#' scored on the same folds. The selected configuration is the one with the
#' lowest held-out perplexity among those that scored on every fold; a
#' configuration with a failed fold is kept in the grid but is never selected.
#' It returns the grid as a \code{transitiontrees_tune} data frame with one row
#' per configuration, lowest held-out perplexity first.
#'
#' The grid compares settings by out-of-sample prediction. To describe
#' the fit of one tree to its own data use \code{\link{model_fit}()} or
#' \code{\link{perplexity}()}; to compare smoothing methods or pruning
#' rules on one data set without cross-validation use
#' \code{\link{compare_smoothing}()} and \code{\link{compare_pruning}()}.
#'
#' @param data Sequences in any format accepted by
#'   \code{\link{context_tree}()}.
#' @param max_depth Maximum depths to try. Default \code{2:5}.
#' @param min_count Minimum context counts to try. Default
#'   \code{c(3, 5, 10)}.
#' @param smoothing Smoothing methods to try: a character vector of method
#'   names, each used with its default settings (default \code{"floor"}),
#'   or a list of specifications such as
#'   \code{list(list("floor", ymin = 0.001), list("floor", ymin = 0.005))}
#'   to vary the settings of one method.
#' @param prune Whether to prune with the \eqn{G^2} criterion; both values
#'   are tried by default (\code{c(FALSE, TRUE)}).
#' @param alpha Significance level of the \eqn{G^2} pruning test. Default
#'   0.05.
#' @param folds Number of cross-validation folds. Default 5.
#' @param seed Seed for the assignment of sequences to folds. Default 1.
#'   The call sets the session's random-number state.
#' @param actor,time,action,order,session,time_threshold Column names and
#'   session gap for a long event log, as in \code{\link{context_tree}()};
#'   the log is reshaped with \code{\link{prepare_input}()} before tuning.
#'   Default \code{NULL} (data already one row per sequence) and 900
#'   seconds.
#'
#' @return An object of class \code{transitiontrees_tune}: a data frame
#'   with one row per configuration, sorted by held-out perplexity (lowest
#'   first), and columns \code{max_depth}, \code{min_count},
#'   \code{smoothing}, \code{prune}, \code{logLik} (held-out
#'   log-likelihood summed over the folds that scored; \code{NA} when none
#'   did), \code{n_scored} (held-out transitions scored),
#'   \code{perplexity}, \code{n_nodes_avg} (mean number of contexts of the
#'   fitted trees) and \code{folds_failed} (folds that could not be fitted
#'   or scored). \code{as.data.frame(x, what = "best")} returns the
#'   selected configuration, \code{summary()} a one-row overview, and
#'   \code{plot()} the perplexity surface.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) fewer sequences than
#'     \code{folds}; \code{smoothing} neither a character vector nor a
#'     list; \code{data} not in a sequence format.}
#'   \item{\code{transitiontrees_fold_failed}}{(warning) at least one
#'     configuration had a fold that could not be fitted or scored; the
#'     message gives the first error.}
#'   \item{\code{transitiontrees_no_best}}{(warning) no configuration
#'     scored on every fold, so none is selected.}
#' }
#'
#' @seealso \code{\link{summary.transitiontrees_tune}},
#'   \code{\link{as.data.frame.transitiontrees_tune}} and
#'   \code{\link{plot.transitiontrees_tune}} for the result;
#'   \code{\link{context_tree}()} and \code{\link{prune_tree}()} to fit the
#'   selected configuration; \code{\link{perplexity}()} and
#'   \code{\link{model_fit}()} for the fit of one tree.
#'
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = c(5, 10),
#'                    prune = c(FALSE, TRUE), folds = 3)
#' tuned
#' summary(tuned)
#' @export
tune_tree <- function(data,
                          max_depth   = 2L:5L,
                          min_count   = c(3L, 5L, 10L),
                          smoothing   = "floor",
                          prune       = c(FALSE, TRUE),
                          alpha       = 0.05,
                          folds       = 5L,
                          seed        = 1L,
                          actor = NULL, time = NULL, action = NULL,
                          order = NULL, session = NULL,
                          time_threshold = 900) {
  ## Long-format data is reshaped first (same as context_tree), so a
  ## long event log is never silently mis-read as one row per sequence.
  data  <- .ct_maybe_reshape(data, actor, time, action, order, session,
                             time_threshold)
  trajs <- .ct_traj(.ct_coerce(data))
  if (length(trajs) < folds)
    .tt_abort("Not enough sequences (", length(trajs), ") for ", folds, " folds.")
  alphabet <- .ct_alphabet(trajs)
  fold_sets <- .pt_make_folds(length(trajs), folds, seed = seed)

  smoothing_grid <- .pt_smoothing_grid(smoothing)
  sm_labels      <- vapply(smoothing_grid, .pt_smoothing_label,
                            character(1))

  base <- expand.grid(max_depth   = as.integer(max_depth),
                      min_count   = as.integer(min_count),
                      smoothing_i = seq_along(smoothing_grid),
                      prune       = as.logical(prune),
                      stringsAsFactors = FALSE)

  scored <- lapply(seq_len(nrow(base)), function(gi) {
    sm <- smoothing_grid[[base$smoothing_i[gi]]]
    cfg <- list(max_depth   = base$max_depth[gi],
                min_count   = base$min_count[gi],
                smoothing   = sm,
                prune       = base$prune[gi],
                alpha       = alpha)
    fold_results <- lapply(fold_sets, function(test_idx) {
      train_idx <- setdiff(seq_along(trajs), test_idx)
      tryCatch(
        .pt_fit_score_fold(train_seqs = trajs[train_idx],
                            test_seqs = trajs[test_idx],
                            alphabet  = alphabet,
                            cfg       = cfg),
        ## A fold that cannot be fitted or scored is recorded, not
        ## dropped: it is counted in `folds_failed`, its message is kept
        ## for the warning below, and the configuration is excluded from
        ## "best".
        error = function(e) list(ll = NA_real_, n = 0L,
                                  n_nodes = NA_integer_,
                                  error = conditionMessage(e))
      )
    })
    failed <- vapply(fold_results,
                     function(fr) is.na(fr$ll) || fr$n == 0L, logical(1))
    folds_failed <- sum(failed)
    ## Only the folds that scored contribute; with none, the held-out
    ## log-likelihood is undefined (NA), not 0.
    ok          <- fold_results[!failed]
    ll_total    <- if (length(ok))
      sum(vapply(ok, `[[`, numeric(1), "ll")) else NA_real_
    n_total     <- if (length(ok))
      sum(vapply(ok, `[[`, numeric(1), "n")) else 0
    n_nodes_avg <- if (length(ok))
      mean(vapply(ok, `[[`, numeric(1), "n_nodes")) else NA_real_
    pp <- if (n_total > 0) exp(-ll_total / n_total) else NA_real_
    first_error <- unlist(lapply(fold_results, `[[`, "error"))
    data.frame(logLik         = ll_total,
               n_scored       = as.integer(n_total),
               perplexity     = pp,
               n_nodes_avg    = n_nodes_avg,
               folds_failed   = as.integer(folds_failed),
               first_error    = if (length(first_error)) first_error[1L]
                                else NA_character_,
               stringsAsFactors = FALSE)
  })
  metrics <- do.call(rbind, scored)

  out <- data.frame(
    max_depth  = base$max_depth,
    min_count  = base$min_count,
    smoothing  = sm_labels[base$smoothing_i],
    prune      = base$prune,
    metrics,
    stringsAsFactors = FALSE
  )
  ord <- order(out$perplexity, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  ## A configuration is only a candidate for "best" if every fold scored
  ## (finite perplexity, no failed folds) — a failed config must never be
  ## reported as the winner. Surface failures rather than masking them.
  n_bad <- sum(out$folds_failed > 0L)
  if (n_bad > 0L) {
    msg <- out$first_error[!is.na(out$first_error)][1L]
    .tt_warn(sprintf(
      "%d of %d tuning configuration(s) had at least one CV fold that failed; their scores use only the folds that succeeded and they are excluded from 'best'.%s",
      n_bad, nrow(out),
      if (is.na(msg)) "" else sprintf(" First error: %s", msg)),
      class = "transitiontrees_fold_failed")
  }
  out$first_error <- NULL
  eligible <- out$folds_failed == 0L & is.finite(out$perplexity)
  best <- if (any(eligible)) out[which(eligible)[1L], , drop = FALSE] else NULL
  if (is.null(best))
    .tt_warn("No tuning configuration scored on every fold; 'best' is NULL.",
             class = "transitiontrees_no_best")
  structure(out, best = best, class = c("transitiontrees_tune", "data.frame"))
}

#' Plot the cross-validated perplexity of a tuning grid
#'
#' @description
#' \code{plot()} draws the held-out perplexity surface of a
#' \code{\link{tune_tree}()} grid: perplexity against maximum depth, one line
#' per minimum count, with one panel per smoothing method (and per pruning
#' choice when both were tried). The selected configuration is marked with a
#' star. The plot shows where added depth stops lowering held-out perplexity and
#' how sensitive that point is to the other settings.
#'
#' @param x A \code{transitiontrees_tune} object from
#'   \code{\link{tune_tree}()}.
#' @param ... Ignored.
#' @return A ggplot object: held-out perplexity (y) against
#'   \code{max_depth} (x), colored by \code{min_count}, faceted by
#'   smoothing and pruning.
#' @seealso \code{\link{tune_tree}()},
#'   \code{\link{summary.transitiontrees_tune}}.
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = c(5, 10),
#'                    prune = FALSE, folds = 3)
#' plot(tuned)
#' @export
plot.transitiontrees_tune <- function(x, ...) {
  best <- attr(x, "best")
  df   <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df$prune_label <- ifelse(df$prune, "pruned", "unpruned")
  df$min_count <- factor(df$min_count)
  ## The star must sit in the facet of the selected configuration only, so
  ## it carries the same facetting column as the grid.
  if (!is.null(best))
    best$prune_label <- ifelse(best$prune, "pruned", "unpruned")

  facet <- if (length(unique(df$prune_label)) > 1L)
    ggplot2::facet_grid(prune_label ~ smoothing) else
    ggplot2::facet_wrap(~ smoothing)

  ggplot2::ggplot(df,
    ggplot2::aes(x = .data$max_depth, y = .data$perplexity,
                 colour = .data$min_count, group = .data$min_count)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_point(data = best,
                        ggplot2::aes(x = .data$max_depth,
                                     y = .data$perplexity),
                        shape = 8, size = 4, colour = "#D55E00",
                        inherit.aes = FALSE) +
    facet +
    ggplot2::labs(
      x = "max_depth", y = "held-out perplexity",
      colour = "min_count",
      title = "context tree CV grid",
      subtitle = "lower is better; star = minimum-perplexity configuration"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Print a tuning grid
#'
#' @description
#' \code{print()} prints the number of configurations, the first \code{n} rows
#' of the grid (lowest held-out perplexity first) and the selected
#' configuration.
#'
#' @param x A \code{transitiontrees_tune} object.
#' @param n Number of configurations to print. Default 10.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{tune_tree}()},
#'   \code{\link{as.data.frame.transitiontrees_tune}} for the full grid.
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:2, min_count = 5,
#'                    prune = FALSE, folds = 3)
#' print(tuned, n = 3)
#' @export
print.transitiontrees_tune <- function(x, n = 10L, ...) {
  cat(sprintf("<transitiontrees_tune>  %d configurations\n", nrow(x)))
  print.data.frame(utils::head(x, n), row.names = FALSE)
  best <- attr(x, "best")
  if (!is.null(best)) {
    cat("\nbest (min perplexity):\n")
    print.data.frame(best, row.names = FALSE)
  }
  invisible(x)
}

#' Convert a tuning grid to a data frame
#'
#' @description
#' \code{as.data.frame()} returns the cross-validation grid, or only the
#' selected configuration, as a plain data frame.
#'
#' @param x A \code{transitiontrees_tune} object.
#' @param row.names,optional Ignored.
#' @param what \code{"grid"} (default) for every configuration, or
#'   \code{"best"} for the selected one.
#' @param ... Ignored.
#' @return A data frame. For \code{what = "grid"}, one row per
#'   configuration, sorted by held-out perplexity, with the columns
#'   described in \code{\link{tune_tree}()}. For \code{what = "best"}, one
#'   row: the configuration with the lowest held-out perplexity among those
#'   that scored on every fold, or zero rows when none did.
#' @seealso \code{\link{tune_tree}()},
#'   \code{\link{summary.transitiontrees_tune}}.
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = 5,
#'                    prune = FALSE, folds = 3)
#' as.data.frame(tuned, what = "best")
#' @export
as.data.frame.transitiontrees_tune <- function(x, row.names = NULL,
                                               optional = FALSE,
                                               what = c("grid", "best"),
                                               ...) {
  what <- match.arg(what)
  grid <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(grid, "best") <- NULL
  if (identical(what, "grid")) return(grid)
  best <- attr(x, "best")
  if (is.null(best)) return(grid[0L, , drop = FALSE])
  best <- as.data.frame(unclass(best), stringsAsFactors = FALSE)
  attr(best, "best") <- NULL
  rownames(best) <- NULL
  best
}

#' Summarize a tuning grid
#'
#' @description
#' \code{summary()} reduces a \code{\link{tune_tree}()} grid to one row: the
#' number of configurations scored, the number with a failed cross-validation
#' fold, and the settings, held-out perplexity and mean size of the selected
#' configuration.
#'
#' @param object A \code{transitiontrees_tune} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.transitiontrees_tune}: a
#'   one-row data frame with columns \code{n_configurations},
#'   \code{n_failed}, \code{max_depth}, \code{min_count},
#'   \code{smoothing}, \code{prune}, \code{perplexity} and
#'   \code{n_nodes_avg}. The last six describe the selected configuration
#'   and are \code{NA} when no configuration scored on every fold.
#' @seealso \code{\link{tune_tree}()},
#'   \code{\link{as.data.frame.transitiontrees_tune}}.
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = 5,
#'                    prune = c(FALSE, TRUE), folds = 3)
#' summary(tuned)
#' @export
summary.transitiontrees_tune <- function(object, ...) {
  .tt_check_class(object, "transitiontrees_tune", "object")
  grid <- as.data.frame(object, what = "grid")
  best <- as.data.frame(object, what = "best")
  has_best <- nrow(best) == 1L
  pick <- function(col, na) if (has_best) best[[col]] else na
  out <- data.frame(
    n_configurations = nrow(grid),
    n_failed         = sum(grid$folds_failed > 0L),
    max_depth        = pick("max_depth", NA_integer_),
    min_count        = pick("min_count", NA_integer_),
    smoothing        = pick("smoothing", NA_character_),
    prune            = pick("prune", NA),
    perplexity       = pick("perplexity", NA_real_),
    n_nodes_avg      = pick("n_nodes_avg", NA_real_),
    stringsAsFactors = FALSE, row.names = NULL)
  class(out) <- c("summary.transitiontrees_tune", "data.frame")
  out
}

#' Print a summary of a tuning grid
#'
#' @description
#' \code{print()} prints the number of configurations and failed configurations,
#' then the selected configuration as a one-row table.
#'
#' @param x A \code{summary.transitiontrees_tune} object.
#' @param digits Significant digits for numeric columns. Default 4.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{summary.transitiontrees_tune}}.
#' @examples
#' tuned <- tune_tree(trajectories, max_depth = 1:2, min_count = 5,
#'                    prune = FALSE, folds = 3)
#' print(summary(tuned), digits = 3)
#' @export
print.summary.transitiontrees_tune <- function(x, digits = 4L, ...) {
  cat(sprintf("<transitiontrees_tune summary>  %d configurations, %d with a failed fold\n",
              x$n_configurations, x$n_failed))
  cat("selected (minimum held-out perplexity):\n")
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df$n_configurations <- NULL
  df$n_failed <- NULL
  print.data.frame(df, row.names = FALSE, digits = digits)
  invisible(x)
}
