# Reference computations shared by the equivalence tests. Nothing here calls transitiontrees: each function
# builds the input both sides receive, runs a reference package (PST, TraMineR, entropy, markovchain, tna), or
# writes a definition out in base R. The package is the other side of each comparison.
#
# PST (Gabadinho & Ritschard 2016) is fitted on a TraMineR state sequence object of the same matrix. PST writes
# a context as "A-B-C" (oldest first) and the root as "e"; transitiontrees writes "A -> B -> C" and "(start)".
# PST's nobs(), predict() and logLik() are S4 methods; transitiontrees' S3 methods of the same generics are
# attached in the test session, so the S4 methods are called through methods::getMethod().

#' Sequences of a first-order Markov chain: the shared synthetic input
#'
#' The transition matrix has rows drawn from U(0.05, 1) and normalised; the first state is uniform. The draws
#' are seeded here, so the input is the same in every run.
#' @param n_seq,n_t Number and length of the sequences.
#' @param k Number of states (named A, B, ...).
#' @param seed The seed of the draws.
#' @return A character matrix, one sequence per row.
ref_chain <- function(n_seq, n_t, k, seed) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  states <- LETTERS[seq_len(k)]
  weights <- matrix(stats::runif(k * k, 0.05, 1), k, k, dimnames = list(states, states))
  transition <- weights / rowSums(weights)
  cumulative <- t(apply(transition, 1L, cumsum))
  first <- sample(states, n_seq, replace = TRUE)
  step <- \(previous, t) states[rowSums(stats::runif(n_seq) > cumulative[previous, , drop = FALSE]) + 1L]
  columns <- Reduce(step, seq_len(n_t - 1L), first, accumulate = TRUE)
  out <- do.call(cbind, columns)
  dimnames(out) <- list(NULL, paste0("T", seq_len(n_t)))
  out
}

#' The configurations the fit is compared on
#'
#' Written out, not drawn: alphabet size 2 to 4, depth 2 to 4, min_count 3 to 8, ymin 0 to 0.005.
ref_configs <- data.frame(
  config = paste0("c", 1:8),
  k = c(2L, 3L, 3L, 4L, 4L, 2L, 3L, 4L),
  n_seq = c(60L, 80L, 120L, 80L, 120L, 50L, 100L, 60L),
  n_t = c(20L, 15L, 25L, 20L, 15L, 25L, 20L, 25L),
  max_depth = c(3L, 3L, 4L, 2L, 3L, 4L, 2L, 3L),
  min_count = c(3L, 5L, 5L, 3L, 8L, 5L, 3L, 5L),
  ymin = c(0.001, 0, 0.005, 0.001, 0, 0.005, 0.001, 0.002),
  seed = c(101L, 202L, 303L, 404L, 505L, 606L, 707L, 808L),
  stringsAsFactors = FALSE)

#' The input of one configuration
#' @param row One row of `ref_configs`.
#' @return A list: `config`, `sequences` (character matrix), `alphabet`, and the fit settings.
ref_input <- function(row) {
  list(config = row$config, sequences = ref_chain(row$n_seq, row$n_t, row$k, row$seed),
       alphabet = LETTERS[seq_len(row$k)], max_depth = row$max_depth, min_count = row$min_count, ymin = row$ymin)
}

#' All configurations' inputs, by name
ref_inputs <- function() {
  inputs <- lapply(seq_len(nrow(ref_configs)), \(i) ref_input(ref_configs[i, ]))
  stats::setNames(inputs, ref_configs$config)
}

#' A TraMineR state sequence object of a matrix (NA = missing; trailing NA are void, as TraMineR's default)
pst_sequences <- function(sequences, alphabet) {
  suppressMessages(TraMineR::seqdef(as.data.frame(sequences, stringsAsFactors = FALSE), alphabet = alphabet))
}

#' PST::pstree() of a matrix
pst_fit <- function(sequences, alphabet, max_depth, min_count, ymin, weights = NULL) {
  states <- if (is.null(weights)) pst_sequences(sequences, alphabet) else
    suppressMessages(TraMineR::seqdef(as.data.frame(sequences, stringsAsFactors = FALSE), alphabet = alphabet,
                                      weights = weights))
  suppressMessages(PST::pstree(states, L = max_depth, nmin = min_count, ymin = ymin))
}

#' PST::pstree() of one input of ref_inputs()
pst_fit_input <- function(input, ymin = input$ymin) {
  pst_fit(input$sequences, input$alphabet, input$max_depth, input$min_count, ymin)
}

#' A PST context ("A-B", "e") in transitiontrees' notation ("A -> B", "(start)")
pst_to_arrow <- function(context) ifelse(context == "e", "(start)", gsub("-", " -> ", context, fixed = TRUE))

#' A transitiontrees context in PST's notation
arrow_to_pst <- function(context) ifelse(context == "(start)", "e", gsub(" -> ", "-", context, fixed = TRUE))

#' The nodes of a PST tree, one row per node (pruned nodes left out)
#'
#' @return A data frame: `pathway` (arrow notation), `depth`, `count` (the node's n).
pst_node_table <- function(pst) {
  levels <- methods::as(pst, "list")
  rows <- lapply(seq_along(levels), \(level) {
    nodes <- Filter(\(node) !any(node@pruned), levels[[level]])
    if (!length(nodes)) return(NULL)
    contexts <- vapply(nodes, \(node) if (level == 1L) "e" else paste(node@path, collapse = "-"), character(1))
    data.frame(pathway = pst_to_arrow(contexts), depth = level - 1L,
               count = vapply(nodes, \(node) sum(node@n), numeric(1)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' The nodes of a PST tree, one row per node and next state: counts and probabilities
#'
#' @return A data frame: `pathway`, `state`, `count`, `prob`.
pst_state_table <- function(pst) {
  levels <- methods::as(pst, "list")
  alphabet <- pst@alphabet
  rows <- lapply(seq_along(levels), \(level) {
    nodes <- Filter(\(node) !any(node@pruned), levels[[level]])
    do.call(rbind, lapply(nodes, \(node) {
      context <- if (level == 1L) "e" else paste(node@path, collapse = "-")
      data.frame(pathway = pst_to_arrow(context), state = alphabet, count = as.numeric(node@counts[1L, ]),
                 prob = as.numeric(node@prob[1L, ]), stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' PST::query() of a context in arrow notation, as a named numeric vector
pst_query <- function(pst, pathway) {
  stats::setNames(as.numeric(suppressMessages(PST::query(pst, arrow_to_pst(pathway)))), pst@alphabet)
}

#' The context PST::query() falls back to (longest suffix in the tree), in arrow notation
pst_matched <- function(pst, pathway) {
  pst_to_arrow(suppressMessages(PST::query(pst, arrow_to_pst(pathway)))@context)
}

#' PST's S4 methods of stats generics (see the note at the top)
pst_logLik <- function(pst) suppressMessages(methods::getMethod("logLik", class(pst)[[1L]])(pst))
pst_nobs <- function(pst) as.numeric(methods::getMethod("nobs", class(pst)[[1L]])(pst))
pst_predict <- function(pst, states, ...) suppressMessages(methods::getMethod("predict", class(pst)[[1L]])(pst, states, ...))
pst_prune <- function(pst, ...) suppressMessages(methods::getMethod("prune", class(pst)[[1L]])(pst, ...))

#' Per-position predicted probabilities of PST::predict(decomp = TRUE), one row per scored position
#'
#' @param pst A PST tree. @param sequences A character matrix (NA = missing).
#' @return A data frame: `sequence_id`, `position`, `predicted_prob`.
pst_positions <- function(pst, sequences) {
  states <- pst_sequences(sequences, pst@alphabet)
  probs <- pst_predict(pst, states, output = "prob", decomp = TRUE)
  out <- data.frame(sequence_id = as.vector(row(probs)), position = as.vector(col(probs)),
                    predicted_prob = as.vector(probs))
  out <- stats::na.omit(out)
  attr(out, "na.action") <- NULL
  rownames(out) <- NULL
  out
}

#' PST's summary() of a tree (an S4 method; base summary() would dispatch first)
pst_summary <- function(pst) methods::getMethod("summary", class(pst)[[1L]])(pst)

#' The free parameters PST reports for a tree
pst_freepar <- function(pst) as.numeric(pst_summary(pst)@freepar)

#' The number of nodes PST reports for a tree (internal nodes + leaves)
pst_n_nodes <- function(pst) as.numeric(pst_summary(pst)@nodes) + as.numeric(pst_summary(pst)@leaves)

#' Kullback-Leibler divergence D(p || q) with entropy::KL.plugin(), in the given unit
ref_kl <- function(p, q, unit = "log2") entropy::KL.plugin(p, q, unit = unit)

#' Shannon entropy with entropy::entropy.plugin(), in the given unit
ref_entropy <- function(p, unit = "log2") entropy::entropy.plugin(p, unit = unit)

#' The immediate suffix of a context ("A -> B -> C" -> "B -> C"; a one-state context -> "(start)")
ref_parent <- function(pathway) {
  parts <- strsplit(pathway, " -> ", fixed = TRUE)
  vapply(parts, \(p) if (length(p) <= 1L) "(start)" else paste(p[-1L], collapse = " -> "), character(1))
}

#' Probabilities of a PST tree as a list keyed by context (arrow notation)
pst_prob_list <- function(pst) {
  states <- pst_state_table(pst)
  lapply(split(states, factor(states$pathway, unique(states$pathway))), \(rows) stats::setNames(rows$prob, rows$state))
}

#' The first `n` rows of a data frame ordered by one column
ref_head_by <- function(data, column, n, decreasing = TRUE) {
  utils::head(data[order(data[[column]], decreasing = decreasing), , drop = FALSE], n)
}

#' Per-context counts and probabilities of a PST tree, as a list keyed by context
#'
#' @return A named list; each entry a list of `counts` and `prob` (vectors in alphabet order).
pst_node_list <- function(pst) {
  states <- pst_state_table(pst)
  groups <- split(states, factor(states$pathway, unique(states$pathway)))
  lapply(groups, \(rows) list(counts = rows$count, prob = rows$prob))
}

#' Leaf pruning written out from its definition, on PST's counts and probabilities
#'
#' A leaf (a context that is no other context's immediate suffix) whose immediate suffix is in the tree is
#' tested against that suffix; the leaves that fail are removed, and the test repeats on the new leaves until
#' none fails. The tests, for a leaf with counts n_s (total n), probabilities p and parent probabilities q:
#' G2: 2 * sum(n_s * log(n_s / (n * q_s))) > cut; KL: KL(p || q) in nats > cut; AIC: 2 (k - 1) - 2 sum(n_s log p_s)
#' < -2 sum(n_s log q_s); BIC: log(n) (k - 1) - 2 sum(n_s log p_s) < -2 sum(n_s log q_s).
#' @param pst A PST tree. @param criterion "G2", "KL", "AIC" or "BIC". @param cut The cutoff of G2 or KL.
#' @return The contexts kept (arrow notation, the root included).
ref_prune <- function(pst, criterion, cut) {
  nodes <- pst_node_list(pst)
  k <- length(pst@alphabet)
  keeps <- function(context) {
    child <- nodes[[context]]
    parent <- nodes[[ref_parent(context)]]
    counts <- child$counts
    seen <- counts > 0
    n <- sum(counts)
    ll_child <- sum(counts[seen] * log(child$prob[seen]))
    ll_parent <- sum(counts[seen] * log(parent$prob[seen]))
    switch(criterion,
           G2 = 2 * sum(counts[seen] * log(counts[seen] / (n * parent$prob[seen]))) > cut,
           KL = ref_kl(child$prob, parent$prob, unit = "log") > cut,
           AIC = 2 * (k - 1) - 2 * ll_child < -2 * ll_parent,
           BIC = log(n) * (k - 1) - 2 * ll_child < -2 * ll_parent)
  }
  step <- function(kept) {
    inner <- setdiff(kept, "(start)")
    leaves <- setdiff(inner, ref_parent(inner))
    leaves <- leaves[ref_parent(leaves) %in% kept]
    dropped <- leaves[!vapply(leaves, keeps, logical(1))]
    if (length(dropped) == 0L) kept else step(setdiff(kept, dropped))
  }
  step(names(nodes))
}

#' A smoother written out from its definition
#'
#' For a context with counts n_s (total n, k states) and parent distribution q (the nearest shorter suffix in
#' the tree, already smoothed; uniform for the root). laplace: (n_s + alpha) / (n + k alpha) (Chen & Goodman 1999,
#' additive). jelinek_mercer: (1 - lambda) n_s / n + lambda q_s (Jelinek & Mercer 1980). witten_bell: the same with
#' lambda = r / (r + n), r the number of distinct states seen (Witten & Bell 1991). kneser_ney: absolute
#' discounting, max(n_s - d, 0) / n + (1 - sum(max(n_s - d, 0)) / n) q_s (Ney, Essen & Kneser 1994; the
#' parent's smoothed distribution, not continuation counts). floor_cap: max(n_s / n, ymin) renormalised.
#' @param counts,parent Count vector and parent distribution. @param method,parameter The smoother and its
#'   hyperparameter (alpha, lambda, discount or ymin).
ref_smoother <- function(counts, parent, method, parameter = NULL) {
  k <- length(counts)
  n <- sum(counts)
  mle <- counts / n
  switch(method,
         laplace = (counts + parameter) / (n + k * parameter),
         jelinek_mercer = (1 - parameter) * mle + parameter * parent,
         witten_bell = {
           r <- sum(counts > 0)
           (1 - r / (r + n)) * mle + r / (r + n) * parent
         },
         kneser_ney = {
           discounted <- pmax(counts - parameter, 0) / n
           discounted + (1 - sum(discounted)) * parent
         },
         floor_cap = pmax(mle, parameter) / sum(pmax(mle, parameter)))
}

#' Smooth every context of a tree, shallow to deep, from PST's (unsmoothed) counts
#' @return A named list of probability vectors, keyed by context.
ref_smooth_tree <- function(pst, method, parameter = NULL) {
  nodes <- pst_node_list(pst)
  k <- length(pst@alphabet)
  contexts <- names(nodes)
  depth <- vapply(contexts, \(context) if (context == "(start)") 0L else length(strsplit(context, " -> ", fixed = TRUE)[[1L]]),
                  integer(1))
  ordered <- contexts[order(depth)]
  Reduce(\(smoothed, context) {
    parent <- if (context == "(start)") rep(1 / k, k) else smoothed[[ref_parent(context)]]
    smoothed[[context]] <- ref_smoother(nodes[[context]]$counts, parent, method, parameter)
    smoothed
  }, ordered, list())
}

#' The probability of every observed state under a table of context distributions
#'
#' Each position is predicted by the longest suffix of its history (at most `max_depth` states) that is a
#' context of the table; the first position by the root. As PST::predict() does for its own tree.
#' @param sequences A character matrix without internal gaps. @param probs A named list of distributions
#'   (arrow-notation keys, "(start)" for the root). @param alphabet The states.
#' @return A numeric vector of probabilities, one per observed state.
ref_position_probs <- function(sequences, probs, alphabet, max_depth) {
  rows <- lapply(seq_len(nrow(sequences)), \(i) as.character(stats::na.omit(sequences[i, ])))
  unlist(lapply(rows, \(states) vapply(seq_along(states), \(t) {
    history <- utils::tail(states[seq_len(t - 1L)], max_depth)
    suffixes <- vapply(seq_along(history), \(start) paste(history[start:length(history)], collapse = " -> "), character(1))
    found <- suffixes[suffixes %in% names(probs)]
    context <- if (length(found)) found[[1L]] else "(start)"
    probs[[context]][[match(states[[t]], alphabet)]]
  }, numeric(1))))
}

#' The rows of a matrix split by a grouping vector, as a named list of matrices
ref_split_rows <- function(sequences, group) {
  lapply(split(seq_len(nrow(sequences)), group), \(rows) sequences[rows, , drop = FALSE])
}

#' A table with constant columns put in front (a table of no rows stays one of no rows)
ref_tag <- function(table, ...) {
  values <- list(...)
  data.frame(lapply(values, rep, length.out = nrow(table)), table, stringsAsFactors = FALSE, check.names = FALSE)
}
