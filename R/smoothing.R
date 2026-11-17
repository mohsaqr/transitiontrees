# ---- Smoothing schemes for transitiontrees ----
#
# Implements:
#   .smooth_floor          - floor MLE at ymin (default "interpolate"
#                            rule; "cap" = original renormalise)
#   .smooth_laplace        - additive-alpha
#   .smooth_kneser_ney     - back-off along suffix path (parent-prob
#                            approximation; see notes)
#   .smooth_witten_bell    - novelty-weighted interpolation
#   .smooth_jelinek_mercer - fixed-lambda interpolation
#   .ct_smooth_dispatch    - switch dispatcher
#   smooth_tree()      - re-smooth a fitted tree without refitting

#' Floor smoothing of a count vector to a probability distribution.
#'
#' Two rules:
#' \itemize{
#'   \item \code{"interpolate"} (default): a distribution
#'     that has at least one zero-count state is shifted toward uniform,
#'     \eqn{p_i = (1 - k\,y_{min}) p_i + y_{min}} (k = alphabet size), so
#'     each zero state lands at exactly \code{ymin}. Distributions with
#'     every state observed are left as the raw MLE.
#'   \item \code{"cap"}: clamp every probability up to \code{ymin} and
#'     renormalise (\code{pmax(p, ymin) / sum(...)}). transitiontrees's original
#'     rule; kept available for back-compatibility.
#' }
#' @noRd
.smooth_floor <- function(counts, parent_prob = NULL, ...,
                          ymin = 0.001, rule = c("interpolate", "cap")) {
  rule <- match.arg(rule)
  k <- length(counts); n <- sum(counts)
  if (n == 0) return(if (is.null(parent_prob)) rep(1 / k, k) else parent_prob)
  p <- counts / n
  if (ymin <= 0) return(p)
  if (rule == "cap") {
    q <- pmax(p, ymin)
    return(q / sum(q))
  }
  ## "interpolate": shift a zero-containing distribution toward uniform.
  ## The coefficient (1 - k*ymin) must stay positive, else probabilities
  ## go negative; require ymin < 1/k.
  if (k * ymin >= 1)
    .tt_abort("floor smoothing with rule = \"interpolate\" needs ymin < 1/k ",
         "(k = alphabet size = ", k, "); got ymin = ", ymin,
         ". Lower ymin, or use rule = \"cap\".")
  if (any(p == 0)) p <- (1 - k * ymin) * p + ymin
  p
}

#' @noRd
.smooth_laplace <- function(counts, parent_prob = NULL, ...,
                            alpha = 1) {
  k <- length(counts)
  (counts + alpha) / (sum(counts) + alpha * k)
}

#' @noRd
.smooth_kneser_ney <- function(counts, parent_prob = NULL, ...,
                               discount = 0.75) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  high   <- pmax(counts - discount, 0) / n
  ## Back-off weight = the probability mass removed by the absolute
  ## discount, computed as the residual so the result sums to exactly
  ## 1. For integer counts this equals discount * n_pos / n (every
  ## positive count is >= 1 >= discount); for fractional (weighted)
  ## counts a cell below `discount` removes only its own mass, which
  ## the residual captures correctly where discount * n_pos / n would
  ## over-count and break the sum-to-1 constraint.
  back_w <- 1 - sum(high)
  high + back_w * parent_prob
}

#' @noRd
.smooth_witten_bell <- function(counts, parent_prob = NULL, ...) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  n_pos  <- sum(counts > 0)
  lambda <- n_pos / (n_pos + n)
  mle    <- counts / n
  (1 - lambda) * mle + lambda * parent_prob
}

#' @noRd
.smooth_jelinek_mercer <- function(counts, parent_prob = NULL, ...,
                                   lambda = 0.5) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  mle <- counts / n
  (1 - lambda) * mle + lambda * parent_prob
}

#' @noRd
.ct_smooth_dispatch <- function(smoothing, counts, parent_prob = NULL) {
  ## `smoothing` is a resolved list from .pt_resolve_smoothing().
  switch(smoothing$method,
    floor          = .smooth_floor(counts, parent_prob,
                                   ymin = smoothing$ymin,
                                   rule = smoothing$rule),
    laplace        = .smooth_laplace(counts, parent_prob,
                                     alpha = smoothing$alpha),
    kneser_ney     = .smooth_kneser_ney(counts, parent_prob,
                                        discount = smoothing$discount),
    witten_bell    = .smooth_witten_bell(counts, parent_prob),
    jelinek_mercer = .smooth_jelinek_mercer(counts, parent_prob,
                                            lambda = smoothing$lambda),
    .tt_abort("Unknown smoothing method: '", smoothing$method, "'."))
}

#' @noRd
.pt_smoothing_defaults <- list(
  floor          = list(method = "floor",          ymin = 0.001,
                        rule = "interpolate"),
  laplace        = list(method = "laplace",        alpha = 1),
  kneser_ney     = list(method = "kneser_ney",     discount = 0.75),
  witten_bell    = list(method = "witten_bell"),
  jelinek_mercer = list(method = "jelinek_mercer", lambda = 0.5)
)

#' @noRd
.pt_resolve_smoothing <- function(spec) {
  ## Accepts a method name (character scalar) or a list whose first
  ## element / "method" entry names the scheme; remaining named entries
  ## override the method's defaults.
  if (is.character(spec) && length(spec) == 1L) {
    if (!spec %in% names(.pt_smoothing_defaults))
      .tt_abort("Unknown smoothing method: '", spec, "'. Choose one of: ",
           paste(names(.pt_smoothing_defaults), collapse = ", "))
    return(.pt_validate_smoothing(.pt_smoothing_defaults[[spec]]))
  }
  if (is.list(spec)) {
    method <- if (!is.null(spec$method)) spec$method else spec[[1L]]
    if (!is.character(method) ||
        !method %in% names(.pt_smoothing_defaults))
      .tt_abort("Smoothing method must be one of: ",
           paste(names(.pt_smoothing_defaults), collapse = ", "))
    nm <- names(spec)
    overrides <- if (is.null(nm)) list() else
      spec[nzchar(nm) & nm != "method"]
    return(.pt_validate_smoothing(
      modifyList(.pt_smoothing_defaults[[method]], overrides)))
  }
  .tt_abort("'smoothing' must be a method name (character) or a list ",
       "(method, ...kwargs).")
}

#' @noRd
.pt_validate_smoothing <- function(sm) {
  scalar_num <- function(x) is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x)
  bad <- function(arg, rule) {
    .tt_abort("Invalid smoothing parameter '", arg, "': expected ",
         rule, ".")
  }
  switch(sm$method,
    floor = {
      if (!scalar_num(sm$ymin) || sm$ymin < 0)
        bad("ymin", "a finite number >= 0")
      if (!is.null(sm$rule) && !(length(sm$rule) == 1L &&
          sm$rule %in% c("interpolate", "cap")))
        bad("rule", 'either "interpolate" or "cap"')
    },
    laplace = {
      if (!scalar_num(sm$alpha) || sm$alpha < 0)
        bad("alpha", "a finite number >= 0")
    },
    kneser_ney = {
      if (!scalar_num(sm$discount) || sm$discount < 0 || sm$discount > 1)
        bad("discount", "a finite number in [0, 1]")
    },
    witten_bell = NULL,
    jelinek_mercer = {
      if (!scalar_num(sm$lambda) || sm$lambda < 0 || sm$lambda > 1)
        bad("lambda", "a finite number in [0, 1]")
    },
    .tt_abort("Unknown smoothing method: '", sm$method, "'."))
  sm
}

#' @noRd
.pt_last_state <- function(ctx) {
  if (identical(ctx, .ROOT)) return(.ROOT_LABEL)
  strsplit(ctx, " -> ", fixed = TRUE)[[1L]][[1L]]
}

#' @noRd
.pt_children_of <- function(tree) {
  if (nrow(tree$edges) == 0L) return(list())
  split(tree$edges$child, tree$edges$parent)
}

#' @noRd
.pt_parent_ctx <- function(ctx) {
  if (identical(ctx, .ROOT)) return(NA_character_)
  parts <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
  if (length(parts) == 1L) return(.ROOT)
  paste(parts[-1L], collapse = " -> ")
}

#' @noRd
.pt_parent_prob <- function(nodes, ctx) {
  ## Walks up the suffix path to the nearest ancestor present in
  ## `nodes`; needed because nmin can drop intermediate ancestors.
  if (identical(ctx, .ROOT)) return(NULL)
  cur <- ctx
  while (TRUE) {
    cur <- .pt_parent_ctx(cur)
    if (is.na(cur)) return(NULL)
    if (cur %in% names(nodes)) return(nodes[[cur]]$prob)
    if (identical(cur, .ROOT)) return(NULL)
  }
}

#' Re-smooth the next-state probabilities of a fitted tree
#'
#' @description
#' \code{smooth_tree()} re-estimates the next-state probabilities of every
#' context of a fitted tree under a new smoothing scheme, without refitting.
#' Smoothing moves part of the probability toward unobserved states or toward
#' the distribution of the shorter context, the parent in the tree (Chen &
#' Goodman, 1999). The next-state distribution of a context is estimated from
#' the continuations observed after it; raw proportions give probability zero to
#' every continuation that was not observed, which makes any new sequence
#' containing such a continuation impossible under the model, and they are
#' unstable in contexts with few occurrences.
#'
#' Counts and the set of contexts stay as they are, and only the probabilities
#' change. Parents are smoothed before their children, so the back-off schemes
#' interpolate with an already smoothed parent. The five schemes are given under
#' Details. It returns a tree of the same class with re-smoothed probabilities.
#'
#' Smoothing changes predictions, perplexity and the \code{divergence}
#' column of \code{\link{tree_pathways}()}, but never which contexts the
#' tree contains; that is set by \code{min_count} and by pruning. To compare
#' schemes on one dataset use \code{\link{compare_smoothing}()} (in-sample)
#' or \code{\link{tune_tree}()} (cross-validated).
#'
#' @details
#' For a context with counts \eqn{n_s} over \eqn{k} states, total \eqn{n}
#' and parent distribution \eqn{q_s} (uniform for the root):
#' \describe{
#'   \item{\code{"floor"}}{With \code{rule = "interpolate"} (default), a
#'     distribution with an unobserved state becomes
#'     \eqn{(1 - k\,y_{min}) p_s + y_{min}}, so that each unobserved state
#'     has probability exactly \code{ymin}; distributions without zeros are
#'     left unchanged. \code{ymin} must be below \eqn{1/k}. With
#'     \code{rule = "cap"}, every probability is raised to at least
#'     \code{ymin} and the vector is renormalized.}
#'   \item{\code{"laplace"}}{Additive smoothing:
#'     \eqn{(n_s + \alpha) / (n + k\alpha)}.}
#'   \item{\code{"kneser_ney"}}{Absolute discounting with back-off to the
#'     parent: \eqn{\max(n_s - d, 0)/n} plus the removed mass times
#'     \eqn{q_s}. The parent's smoothed distribution is used for the
#'     back-off, not the continuation counts of the original method
#'     (Kneser & Ney, 1995).}
#'   \item{\code{"witten_bell"}}{Interpolation with the parent,
#'     \eqn{(1 - \lambda) n_s / n + \lambda q_s}, with
#'     \eqn{\lambda = r / (r + n)} and \eqn{r} the number of distinct
#'     states observed after the context (Witten & Bell, 1991).}
#'   \item{\code{"jelinek_mercer"}}{Interpolation with the parent at a fixed
#'     weight: \eqn{(1 - \lambda) n_s / n + \lambda q_s}.}
#' }
#'
#' @param tree A \code{transitiontrees} object, or a
#'   \code{transitiontrees_group}, whose trees are each re-smoothed.
#' @param smoothing A method name, or a list of the method name and its
#'   hyperparameters, such as \code{list("kneser_ney", discount = 0.5)}.
#'   Methods and default hyperparameters: \code{"floor"} (default;
#'   \code{ymin = 0.001}, \code{rule = "interpolate"}), \code{"laplace"}
#'   (\code{alpha = 1}), \code{"kneser_ney"} (\code{discount = 0.75}),
#'   \code{"witten_bell"} (none) and \code{"jelinek_mercer"}
#'   (\code{lambda = 0.5}).
#'
#' @return An object of the same class as \code{tree} with re-smoothed
#'   probabilities; counts and contexts are unchanged.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree,
#' the method is unknown, a hyperparameter is outside its range, or
#' \code{ymin} is not below \eqn{1/k} under \code{rule = "interpolate"}.
#'
#' @references
#' Begleiter, R., El-Yaniv, R., & Yona, G. (2004). On prediction using
#' variable order Markov models. \emph{Journal of Artificial Intelligence
#' Research}, 22, 385-421. \doi{10.1613/jair.1491}
#'
#' Chen, S. F., & Goodman, J. (1999). An empirical study of smoothing
#' techniques for language modeling. \emph{Computer Speech and Language},
#' 13(4), 359-393. \doi{10.1006/csla.1999.0128}
#'
#' Kneser, R., & Ney, H. (1995). Improved backing-off for M-gram language
#' modeling. In \emph{1995 International Conference on Acoustics, Speech,
#' and Signal Processing} (Vol. 1, pp. 181-184).
#' \doi{10.1109/ICASSP.1995.479394}
#'
#' Witten, I. H., & Bell, T. C. (1991). The zero-frequency problem:
#' Estimating the probabilities of novel events in adaptive text
#' compression. \emph{IEEE Transactions on Information Theory}, 37(4),
#' 1085-1094. \doi{10.1109/18.87000}
#'
#' @seealso \code{\link{compare_smoothing}()} to compare schemes;
#'   \code{\link{context_tree}()}, whose \code{smoothing} argument applies a
#'   scheme at fit time; \code{\link{tune_tree}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' smooth_tree(tree, "kneser_ney")
#' smooth_tree(tree, list("laplace", alpha = 0.5))
#' @export
smooth_tree <- function(tree, smoothing = "floor") {
  ## A transitiontrees_group re-smooths each member, preserving the wrapper.
  if (inherits(tree, "transitiontrees_group")) {
    out <- lapply(tree, smooth_tree, smoothing = smoothing)
    return(structure(out, class = class(tree),
                     group = attr(tree, "group")))
  }
  .tt_check_class(tree, "transitiontrees", "tree")
  sm <- .pt_resolve_smoothing(smoothing)

  ## Re-smooth shallow-to-deep: this loop carries a true sequential
  ## dependency (each node reads its already re-smoothed parent from
  ## new_nodes), so it cannot be vectorised into a single vapply.
  ord <- order(vapply(tree$nodes, function(x) x$depth, integer(1)))
  new_nodes <- tree$nodes
  for (ctx in names(tree$nodes)[ord]) {
    parent_prob <- .pt_parent_prob(new_nodes, ctx)
    new_nodes[[ctx]]$prob <- .ct_smooth_dispatch(
      smoothing   = sm,
      counts      = new_nodes[[ctx]]$counts,
      parent_prob = parent_prob
    )
  }
  tree$nodes     <- new_nodes
  tree$smoothing <- sm
  tree
}

#' Compare smoothing schemes on one dataset
#'
#' @description
#' \code{compare_smoothing()} fits a tree under each of several smoothing
#' schemes, with every other setting held fixed, and reports one row per scheme
#' with the size of the tree and its in-sample perplexity. Given an already
#' fitted tree, it re-smooths that tree under each scheme instead, which keeps
#' its contexts, for example after pruning.
#'
#' The number of contexts is the same for every scheme, because smoothing
#' changes probabilities, not which contexts exist. In-sample perplexity
#' favors schemes that stay close to the observed proportions and therefore
#' is not a criterion for choosing one; \code{\link{tune_tree}()} compares
#' schemes on held-out data. See \code{\link{smooth_tree}()} for the
#' definitions of the schemes.
#'
#' @param data Sequences in any form accepted by \code{\link{context_tree}()},
#'   fitted afresh under each scheme, or a fitted \code{transitiontrees},
#'   re-smoothed under each scheme.
#' @param smoothing Names of the schemes to compare. Default all five:
#'   \code{"floor"}, \code{"laplace"}, \code{"kneser_ney"},
#'   \code{"witten_bell"}, \code{"jelinek_mercer"}, each with its default
#'   hyperparameters.
#' @param ... Further arguments to \code{\link{context_tree}()}, such as
#'   \code{max_depth} and \code{min_count}, held fixed across schemes. Not
#'   used when \code{data} is a fitted tree.
#'
#' @return A \code{data.frame} with one row per scheme, in the order of
#'   \code{smoothing}, and columns \code{smoothing}, \code{n_nodes} (number
#'   of contexts) and \code{perplexity} (in-sample).
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{smoothing} is not a non-empty
#' character vector of known methods, or when \code{data} cannot be fitted
#' (see \code{\link{context_tree}()}).
#'
#' @references
#' Chen, S. F., & Goodman, J. (1999). An empirical study of smoothing
#' techniques for language modeling. \emph{Computer Speech and Language},
#' 13(4), 359-393. \doi{10.1006/csla.1999.0128}
#'
#' @seealso \code{\link{smooth_tree}()} for the schemes;
#'   \code{\link{tune_tree}()} for cross-validated selection;
#'   \code{\link{compare_pruning}()} for the analogous comparison of pruning
#'   criteria.
#'
#' @examples
#' compare_smoothing(trajectories, max_depth = 3)
#' compare_smoothing(trajectories, smoothing = c("floor", "kneser_ney"),
#'                   max_depth = 2, min_count = 10)
#' @export
compare_smoothing <- function(data,
                              smoothing = c("floor", "laplace", "kneser_ney",
                                            "witten_bell", "jelinek_mercer"),
                              ...) {
  if (!is.character(smoothing) || length(smoothing) < 1L)
    .tt_abort("'smoothing' must be a non-empty character vector of method names.")
  ## A fitted tree is re-smoothed (topology frozen, no re-count); raw
  ## data is fitted afresh under each scheme.
  fits <- if (inherits(data, "transitiontrees"))
    lapply(smoothing, function(s) smooth_tree(data, s))
  else
    lapply(smoothing, function(s) context_tree(data, smoothing = s, ...))
  data.frame(
    smoothing  = smoothing,
    n_nodes    = vapply(fits, n_nodes, integer(1)),
    perplexity = vapply(fits, perplexity, numeric(1)),
    stringsAsFactors = FALSE,
    row.names  = NULL)
}
