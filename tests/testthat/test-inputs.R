# Tests for per-sequence weights and TraMineR-style missing codes

mk_mat <- function(n = 30, len = 12, seed = 1,
                   states = c("A","B","C")) {
  set.seed(seed)
  matrix(sample(states, n * len, replace = TRUE), n, len)
}

test_that("weights = rep(1, n) gives byte-identical fit to unweighted", {
  m <- mk_mat()
  tr1 <- context_tree(m, max_depth = 2L, min_count = 3L)
  tr2 <- context_tree(m, max_depth = 2L, min_count = 3L,
                      weights = rep(1, nrow(m)))
  ## Unweighted is integer; weighted is numeric — compare values.
  expect_equal(names(tr1$nodes), names(tr2$nodes))
  for (ctx in names(tr1$nodes)) {
    expect_equal(as.numeric(tr1$nodes[[ctx]]$counts),
                 as.numeric(tr2$nodes[[ctx]]$counts))
    expect_equal(tr1$nodes[[ctx]]$prob,
                 tr2$nodes[[ctx]]$prob)
  }
})

test_that("weighted fit equals unweighted fit on row-duplicated data", {
  ## Duplicating a row twice should be equivalent to weighting that row 3.
  set.seed(11)
  m <- mk_mat()
  dup_idx <- c(seq_len(nrow(m)), 1L, 1L)   # row 1 appears 3x in total
  m_dup <- m[dup_idx, , drop = FALSE]
  tr_dup <- context_tree(m_dup, max_depth = 2L, min_count = 1L)

  w <- rep(1, nrow(m)); w[1L] <- 3
  tr_w   <- context_tree(m, max_depth = 2L, min_count = 1L, weights = w)

  expect_equal(sort(names(tr_dup$nodes)), sort(names(tr_w$nodes)))
  for (ctx in names(tr_dup$nodes)) {
    expect_equal(as.numeric(tr_dup$nodes[[ctx]]$counts),
                 as.numeric(tr_w$nodes[[ctx]]$counts))
  }
})

test_that("weights validation errors on mismatched length", {
  m <- mk_mat()
  expect_error(context_tree(m, max_depth = 1L, min_count = 1L,
                            weights = c(1, 2, 3)),
               "length equal to number of input")
})

test_that("weights validation errors on negative values", {
  m <- mk_mat()
  w <- rep(1, nrow(m)); w[1L] <- -1
  expect_error(context_tree(m, max_depth = 1L, min_count = 1L, weights = w),
               "non-negative")
})

test_that("TraMineR void '%' and missing '*' are treated as gaps, not states", {
  ## A matrix using TraMineR's conventional codes must yield the same
  ## alphabet and fit as the same data with those cells blanked to NA.
  m_codes <- rbind(
    c("A", "B", "%", "*"),
    c("B", "A", "A", "%"),
    c("A", "A", "B", "B"),
    c("B", "%", "*", "%")
  )
  m_blank <- m_codes
  m_blank[m_blank %in% c("%", "*")] <- NA

  tr_codes <- context_tree(m_codes, max_depth = 1L, min_count = 1L)
  tr_blank <- context_tree(m_blank, max_depth = 1L, min_count = 1L)

  expect_equal(tr_codes$alphabet, c("A", "B"))   # no "%" / "*" state
  expect_equal(names(tr_codes$nodes), names(tr_blank$nodes))
  for (ctx in names(tr_codes$nodes))
    expect_equal(as.numeric(tr_codes$nodes[[ctx]]$counts),
                 as.numeric(tr_blank$nodes[[ctx]]$counts))
})

test_that("'%' and '*' are dropped from list input too", {
  seqs <- list(c("A", "%", "B", "A"), c("B", "*", "A"), c("A", "B", "B"))
  tr <- context_tree(seqs, max_depth = 1L, min_count = 1L)
  expect_equal(tr$alphabet, c("A", "B"))
})
