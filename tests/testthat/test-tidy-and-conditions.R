# ---- Tidy prediction/simulation output, classed conditions, tune fixes ----

fit_fixture <- function() context_tree(trajectories, max_depth = 2L)

test_that("predict(tidy = TRUE) carries the same numbers as the matrix", {
  tr <- fit_fixture()
  hist <- list(c("Active", "Average"), "Disengaged", character(0))
  mat <- predict(tr, newdata = hist)
  cls <- predict(tr, newdata = hist, type = "class")
  td  <- predict(tr, newdata = hist, tidy = TRUE)
  expect_s3_class(td, "data.frame")
  expect_named(td, c("id", "history", "context", "predicted", "probability",
                     colnames(mat)))
  expect_identical(nrow(td), 3L)
  expect_identical(td$predicted, cls)
  expect_equal(as.matrix(td[, colnames(mat)]), mat,
               ignore_attr = TRUE, tolerance = 0)
  expect_equal(td$probability, apply(mat, 1L, max))
  expect_identical(td$history,
                   c("Active -> Average", "Disengaged", "(start)"))
  ## each row's distribution sums to one (invariant)
  expect_equal(rowSums(td[, colnames(mat)]), rep(1, 3L),
               tolerance = sqrt(.Machine$double.eps))
  ## default output is unchanged
  expect_true(is.matrix(predict(tr, newdata = hist)))
  expect_error(predict(tr, newdata = hist, tidy = NA))
})

test_that("generate_sequences(tidy = TRUE) is the matrix as a wide data.frame", {
  tr <- fit_fixture()
  set.seed(7); m <- generate_sequences(tr, n = 4L, length = 6L)
  set.seed(7); d <- generate_sequences(tr, n = 4L, length = 6L, tidy = TRUE)
  expect_s3_class(d, "data.frame")
  expect_named(d, paste0("T", 1:6))
  expect_identical(unname(as.matrix(d)), unname(m))
  s <- simulate(tr, nsim = 3L, seed = 3L, length = 5L, tidy = TRUE)
  expect_s3_class(s, "data.frame")
  expect_identical(dim(s), c(3L, 5L))
  ## a simulated frame can be refitted directly
  expect_s3_class(context_tree(s, max_depth = 1L), "transitiontrees")
})

test_that("input errors are classed", {
  tr <- fit_fixture()
  expect_error(context_tree(trajectories, max_depth = -1),
               class = "transitiontrees_bad_input")
  expect_error(context_tree(trajectories, max_depth = -1),
               class = "transitiontrees_error")
  expect_error(tree_pathways(list()), class = "transitiontrees_bad_input")
  expect_error(generate_sequences(tr, n = 0), class = "transitiontrees_bad_input")
  expect_error(query_pathway(tr, "Active", next_state = "Nope"),
               class = "transitiontrees_bad_input")
  other <- context_tree(list(c("x", "y", "x")), max_depth = 1L)
  expect_error(tree_distance(tr, other),
               class = "transitiontrees_alphabet_mismatch")
  w <- context_tree(trajectories, max_depth = 1L,
                    weights = rep(c(1, 2), 69))
  expect_error(bootstrap_pathways(w, iter = 5L),
               class = "transitiontrees_weighted_unsupported")
  expect_error(plot_pathways(tr, min_count = 1e6),
               class = "transitiontrees_empty_result")
})

test_that("tune_tree() names the grid column min_count", {
  set.seed(1)
  m  <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tg <- tune_tree(m, max_depth = 1L, min_count = c(2L, 5L), prune = FALSE,
                  folds = 3L)
  expect_true("min_count" %in% names(tg))
  expect_false("nmin" %in% names(tg))
  expect_setequal(tg$min_count, c(2L, 5L))
  expect_s3_class(plot(tg), "ggplot")
})

test_that("a configuration whose folds all fail has NA logLik, not 0", {
  set.seed(1)
  m  <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  ## ymin = 0.5 >= 1/3 is invalid for the interpolating floor: every fold
  ## of that configuration errors
  expect_warning(
    tg <- tune_tree(m, max_depth = 1L, min_count = 2L, prune = FALSE,
                    smoothing = list(list("floor", ymin = 0.001),
                                     list("floor", ymin = 0.5)),
                    folds = 3L),
    class = "transitiontrees_fold_failed")
  bad <- subset(as.data.frame(tg), folds_failed > 0L)
  expect_identical(nrow(bad), 1L)
  expect_identical(bad$folds_failed, 3L)
  expect_true(is.na(bad$logLik))
  expect_true(is.na(bad$perplexity))
  expect_identical(bad$n_scored, 0L)
  ## the failed configuration is never selected
  expect_identical(as.data.frame(tg, what = "best")$folds_failed, 0L)
  expect_identical(summary(tg)$n_failed, 1L)
})

test_that("no configuration scoring everywhere warns with a class", {
  set.seed(1)
  m  <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  expect_warning(
    expect_warning(
      tg <- tune_tree(m, max_depth = 1L, min_count = 2L, prune = FALSE,
                      smoothing = list(list("floor", ymin = 0.5)),
                      folds = 3L),
      class = "transitiontrees_no_best"),
    class = "transitiontrees_fold_failed")
  expect_identical(nrow(as.data.frame(tg, what = "best")), 0L)
  expect_true(is.na(summary(tg)$perplexity))
})

test_that("forward-trajectory probabilities come from query_pathway()", {
  tr <- fit_fixture()
  td <- transitiontrees:::.trajectory_data(tr, min_count = 4L)
  expect_false(anyNA(td$eprob[td$node != "(start)"]))
  one <- td[td$depth == 2L, , drop = FALSE][1L, ]
  parent_last <- sub(".* -> ", "", one$parent)
  expect_equal(one$eprob,
               unname(query_pathway(tr, parent_last)[[one$last]]))
})

test_that("plot_pruning() labels are ASCII", {
  tr <- context_tree(trajectories, max_depth = 2L)
  p <- plot_pruning(tr, "Active -> Average")
  labs <- c(p$labels$title, p$labels$subtitle, levels(p$data$panel))
  expect_false(any(grepl("[^ -~\n]", labs)))
})

test_that("a tna_data object from tna::prepare_data() fits the same tree as the log", {
  skip_if_not_installed("tna")
  td <- tna::prepare_data(group_regulation_long, actor = "Actor",
                          time = "Time", action = "Action")
  a <- context_tree(td, max_depth = 2L)
  b <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                    action = "Action", max_depth = 2L)
  expect_identical(sort(a$alphabet), sort(b$alphabet))
  expect_equal(model_fit(a), model_fit(b))
  expect_equal(tree_distance(a, b), 0)
})

test_that("the tuning plot marks the selected configuration in its own facet", {
  set.seed(1)
  m  <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tg <- tune_tree(m, max_depth = 1:2, min_count = 2L, prune = c(FALSE, TRUE),
                  folds = 3L)
  p <- plot(tg)
  star <- ggplot2::layer_data(p, 3L)
  expect_identical(nrow(star), 1L)
})

test_that("common_pathways(depth =) returns clean row names", {
  tr <- context_tree(trajectories, max_depth = 2L)
  out <- common_pathways(tr, depth = 2L, top = 3L)
  expect_identical(rownames(out), as.character(seq_len(nrow(out))))
  expect_true(all(out$depth == 2L))
})

test_that("predict() reads arrow-notation histories and warns on unknown states", {
  tr <- context_tree(trajectories, max_depth = 2L)
  a <- predict(tr, newdata = list("Active -> Average"))
  b <- predict(tr, newdata = list(c("Active", "Average")))
  expect_identical(a, b)
  expect_identical(predict(tr, newdata = "Active -> Average"),
                   predict(tr, newdata = c("Active", "Average")))
  expect_warning(predict(tr, newdata = list(c("Active", "Nope"))),
                 class = "transitiontrees_unknown_state")
})
