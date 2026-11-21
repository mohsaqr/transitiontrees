# tune_tree(): cross-validated choice of depth and min_count. With as many folds as sequences (leave one out),
# every fold holds one sequence whatever the random permutation, so the result does not depend on the draws:
# for each configuration, PST::pstree() fitted without sequence i scores sequence i (PST::predict(decomp = TRUE)),
# and the held-out log-likelihood, its count and the perplexity exp(-logLik / n) are summed over i. With fewer
# folds the split is drawn inside the function with no argument to supply it, and is not compared.

test_that("tune_tree(folds = number of sequences) equals leave-one-out PST fits", {
  skip_without_reference(c("PST", "TraMineR"))
  sequences <- ref_chain(30L, 15L, 3L, 61L)
  alphabet <- LETTERS[1:3]
  tuned <- tune_tree(sequences, max_depth = 2:3, min_count = 3, smoothing = "floor", prune = FALSE,
                     folds = nrow(sequences), seed = 1)
  grid <- subset(as.data.frame(tuned, what = "grid"), select = c(max_depth, min_count, logLik, n_scored, perplexity, n_nodes_avg))
  theirs <- do.call(rbind, lapply(2:3, \(depth) {
    folds <- lapply(seq_len(nrow(sequences)), \(i) {
      pst <- pst_fit(sequences[-i, , drop = FALSE], alphabet, max_depth = depth, min_count = 3, ymin = 0.001)
      positions <- pst_positions(pst, sequences[i, , drop = FALSE])
      c(logLik = sum(log(positions$predicted_prob)), n = nrow(positions), nodes = pst_n_nodes(pst))
    })
    totals <- Reduce(`+`, folds)
    data.frame(max_depth = depth, min_count = 3, logLik = totals[["logLik"]], n_scored = totals[["n"]],
               perplexity = exp(-totals[["logLik"]] / totals[["n"]]), n_nodes_avg = totals[["nodes"]] / length(folds))
  }))
  expect_reference(grid, theirs, reference = "PST", covers = c("tune_tree()", "as.data.frame.transitiontrees_tune()"),
                   quantity = "leave-one-out held-out logLik, count, perplexity and mean nodes", tolerance = 1e-10,
                   measure = "relative", key = c("max_depth", "min_count"))
  # The selected configuration: the lowest leave-one-out perplexity.
  expect_reference(subset(as.data.frame(summary(tuned)), select = c(max_depth, min_count, perplexity, n_nodes_avg)),
                   subset(ref_head_by(theirs, "perplexity", n = 1L, decreasing = FALSE),
                          select = c(max_depth, min_count, perplexity, n_nodes_avg)),
                   reference = "PST", covers = "summary.transitiontrees_tune()",
                   quantity = "selected configuration and its perplexity", tolerance = 1e-10, measure = "relative",
                   key = c("max_depth", "min_count"))
})
