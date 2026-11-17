# ---- Classed conditions ----
#
# Every error a user can reasonably catch is raised with a class, so that
# callers and tests can test the condition class instead of message text.
# Each error carries its specific class, then `transitiontrees_error`;
# each warning its specific class, then `transitiontrees_warning`.
#
#   transitiontrees_bad_input             an argument or input of the wrong
#                                          type, shape or value (the default)
#   transitiontrees_alphabet_mismatch     two trees with different alphabets
#   transitiontrees_no_data               a tree without its sequences, where
#                                          resampling needs them
#   transitiontrees_weighted_unsupported  a weighted tree passed to a
#                                          resampling or permutation verb
#   transitiontrees_empty_result          nothing left to plot or report
#   transitiontrees_missing_package       a Suggests package that is needed
#                                          but not installed
#   transitiontrees_fold_failed (warning) a tuning configuration with at least
#                                          one cross-validation fold that failed
#   transitiontrees_no_best (warning)     no tuning configuration scored on
#                                          every fold
#   transitiontrees_unknown_state (warning) a prediction history with a state
#                                          outside the alphabet
#   transitiontrees_degenerate_permutation (warning) a stratified permutation
#                                          whose labels cannot move

#' @noRd
.tt_abort <- function(..., class = "transitiontrees_bad_input") {
  stop(errorCondition(paste0(...),
                      class = c(class, "transitiontrees_error"),
                      call = NULL))
}

#' @noRd
.tt_warn <- function(..., class) {
  warning(warningCondition(paste0(...),
                           class = c(class, "transitiontrees_warning"),
                           call = NULL))
}

#' Check that an argument inherits from a package class
#' @noRd
.tt_check_class <- function(x, cls, arg) {
  if (!inherits(x, cls))
    .tt_abort(sprintf("`%s` must be a %s object; got an object of class <%s>.",
                      arg, cls, paste(class(x), collapse = "/")))
  invisible(TRUE)
}
