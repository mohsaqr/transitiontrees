# ---- Calibration against hand-computed values ----
#
# Two sequences, A B A B and A A B, fitted at depth 1 with no smoothing
# (floor ymin = 0), so every count and probability can be worked out by
# hand:
#   root: every position counted once      A = 4, B = 3   -> 4/7, 3/7
#   "A" : next state after each A          A = 1, B = 3   -> 1/4, 3/4
#   "B" : next state after each non-final B A = 1, B = 0  -> 1, 0
# The log-likelihood scores each position with the longest matching
# context: the first position of each sequence with the root, every later
# position with the context of its predecessor:
#   seq 1: P(A) P(B|A) P(A|B) P(B|A) = 4/7 * 3/4 * 1 * 3/4
#   seq 2: P(A) P(A|A) P(B|A)        = 4/7 * 1/4 * 3/4

hand_seqs <- list(c("A", "B", "A", "B"), c("A", "A", "B"))
hand_tree <- function() {
  context_tree(hand_seqs, max_depth = 1L, min_count = 1L,
               smoothing = list("floor", ymin = 0))
}

test_that("context counts equal the hand counts", {
  tr <- hand_tree()
  expect_identical(sort(names(tr$nodes)), sort(c("<root>", "A", "B")))
  expect_identical(as.integer(tr$nodes[["<root>"]]$counts), c(4L, 3L))
  expect_identical(as.integer(tr$nodes[["A"]]$counts), c(1L, 3L))
  expect_identical(as.integer(tr$nodes[["B"]]$counts), c(1L, 0L))
})

test_that("unsmoothed probabilities are the maximum-likelihood estimates", {
  tr <- hand_tree()
  expect_equal(query_pathway(tr, NULL), c(A = 4 / 7, B = 3 / 7),
               tolerance = sqrt(.Machine$double.eps))
  expect_equal(query_pathway(tr, "A"), c(A = 1 / 4, B = 3 / 4),
               tolerance = sqrt(.Machine$double.eps))
  expect_equal(query_pathway(tr, "B"), c(A = 1, B = 0),
               tolerance = sqrt(.Machine$double.eps))
})

test_that("logLik and nobs equal the hand values", {
  tr <- hand_tree()
  hand_ll <- log(4 / 7) + log(3 / 4) + log(1) + log(3 / 4) +
    log(4 / 7) + log(1 / 4) + log(3 / 4)
  expect_equal(as.numeric(logLik(tr)), hand_ll,
               tolerance = sqrt(.Machine$double.eps))
  expect_identical(as.integer(nobs(tr)), 7L)
  expect_equal(perplexity(tr), exp(-hand_ll / 7),
               tolerance = sqrt(.Machine$double.eps))
})

test_that("state_frequencies() returns the hand counts and proportions", {
  sf <- state_frequencies(hand_tree())
  expect_s3_class(sf, "data.frame")
  expect_named(sf, c("state", "count", "proportion"))
  expect_identical(sf$state, c("A", "B"))
  expect_equal(sf$count, c(4, 3))
  expect_equal(sf$proportion, c(4 / 7, 3 / 7),
               tolerance = sqrt(.Machine$double.eps))
})

test_that("state_frequencies() matches table() on bundled long data", {
  tr <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
                     action = "Action", max_depth = 1L)
  sf <- state_frequencies(tr)
  tab <- table(group_regulation_long$Action)
  expect_equal(sf$count, as.numeric(tab[sf$state]))
  expect_equal(sum(sf$proportion), 1, tolerance = sqrt(.Machine$double.eps))
  ## sorted by count, decreasing
  expect_false(is.unsorted(rev(sf$count)))
})

test_that("state_frequencies() on a group has one block per group", {
  grp <- context_tree(trajectories, group = rep(c("early", "late"), each = 69),
                      max_depth = 1L)
  sf <- state_frequencies(grp)
  expect_named(sf, c("group", "state", "count", "proportion"))
  expect_identical(unique(sf$group), c("early", "late"))
  sums <- tapply(sf$proportion, sf$group, sum)
  expect_equal(as.numeric(sums), c(1, 1), tolerance = sqrt(.Machine$double.eps))
  expect_error(state_frequencies(list(1)),
               class = "transitiontrees_bad_input")
})
