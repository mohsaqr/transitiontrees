# ---- summary() objects and as.data.frame(what = ) accessors ----

grp_fixture <- function() {
  context_tree(trajectories, group = rep(c("early", "late"), each = 69),
               max_depth = 2L)
}

test_that("summary() of a transitiontrees_group is one row per group", {
  grp <- grp_fixture()
  s <- summary(grp)
  expect_s3_class(s, "summary.transitiontrees_group")
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 2L)
  expect_identical(s$group, c("early", "late"))
  expect_named(s, c("group", "n_seq", "n_obs", "n_nodes", "max_depth",
                    "pruned", "logLik", "df", "nobs", "AIC", "BIC",
                    "perplexity"))
  ## agrees with the per-group accessors
  expect_identical(s$n_nodes, unname(n_nodes(grp)))
  mf <- model_fit(grp)
  expect_equal(s$logLik, mf$logLik)
  expect_equal(s$perplexity, mf$perplexity)
  expect_output(print(s), "transitiontrees_group summary")
  expect_invisible(print(s))
})

test_that("summary() of a two-tree comparison is one row and prints", {
  cmp <- compare_trees(grp_fixture(), iter = 20L)
  s <- summary(cmp)
  expect_s3_class(s, "summary.transitiontrees_comparison")
  expect_identical(nrow(s), 1L)
  expect_named(s, c("distance", "null_mean", "null_sd", "null_q95",
                    "p_value", "iter", "n_pathways"))
  nul <- as.data.frame(cmp, what = "null")
  expect_named(nul, c("permutation", "distance"))
  expect_identical(nrow(nul), 20L)
  expect_equal(s$null_mean, mean(nul$distance))
  expect_equal(s$distance, tree_distance(grp_fixture()[["early"]],
                                         grp_fixture()[["late"]]))
  expect_identical(as.data.frame(cmp), as.data.frame(cmp, what = "pathways"))
  expect_output(print(s), "comparison summary")
})

test_that("summary() of a group comparison returns an object, prints nothing", {
  set.seed(1)
  cmp <- compare_groups(grp_fixture(), iter = 49L)
  expect_silent(s <- summary(cmp))
  expect_s3_class(s, "summary.transitiontrees_group_comparison")
  expect_identical(nrow(s), 1L)
  expect_identical(s$groups, "early vs late")
  om <- as.data.frame(cmp, what = "omnibus")
  expect_named(om, c("axis", "statistic", "value", "p_value"))
  expect_equal(s$behavioral_p, om$p_value[om$axis == "behavioral"])
  expect_equal(s$usage_p, om$p_value[om$axis == "usage"])
  pw <- as.data.frame(cmp, what = "pathways")
  expect_identical(s$n_flips, sum(pw$flips))
  expect_identical(s$n_behavioral, sum(pw$jsd_padj < 0.05, na.rm = TRUE))
  ## a looser alpha can only count more pathways (monotone in alpha)
  expect_gte(summary(cmp, alpha = 0.5)$n_usage, s$n_usage)
  expect_output(print(s), "behavioral omnibus")
  expect_error(summary(cmp, alpha = 2))
})

test_that("as.data.frame(what = 'distance') is one row per pair of groups", {
  set.seed(2)
  three <- context_tree(trajectories,
                        group = rep(c("a", "b", "c"), each = 46),
                        max_depth = 1L)
  cmp <- compare_groups(three, iter = 19L)
  d <- as.data.frame(cmp, what = "distance")
  expect_named(d, c("group_a", "group_b", "distance"))
  expect_identical(nrow(d), 3L)
  expect_equal(d$distance[d$group_a == "a" & d$group_b == "b"],
               tree_distance(three[["a"]], three[["b"]]))
  expect_true(all(d$distance >= 0))
})

test_that("summary() and as.data.frame() of a tuning grid", {
  set.seed(1)
  m  <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tg <- tune_tree(m, max_depth = 1:2, min_count = c(2L, 5L), prune = FALSE,
                  folds = 3L)
  grid <- as.data.frame(tg)
  expect_false(inherits(grid, "transitiontrees_tune"))
  expect_null(attr(grid, "best"))
  expect_identical(nrow(grid), 4L)
  best <- as.data.frame(tg, what = "best")
  expect_identical(nrow(best), 1L)
  expect_equal(best$perplexity, min(grid$perplexity))
  s <- summary(tg)
  expect_s3_class(s, "summary.transitiontrees_tune")
  expect_identical(s$n_configurations, 4L)
  expect_identical(s$n_failed, 0L)
  expect_identical(s$max_depth, best$max_depth)
  expect_identical(s$min_count, best$min_count)
  expect_output(print(s), "selected")
})

test_that("as.data.frame(boot, what = 'resamples') is one row per resample x pathway", {
  set.seed(1)
  seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  boot <- bootstrap_pathways(context_tree(seqs, max_depth = 1L), iter = 10L)
  sm <- as.data.frame(boot)
  expect_identical(sm, summary(boot))
  rs <- as.data.frame(boot, what = "resamples", stat = "G2")
  expect_named(rs, c("resample", "pathway", "stat", "value"))
  expect_identical(nrow(rs), 10L * ncol(boot$M_G2))
  expect_identical(unique(rs$stat), "G2")
  ## values are the stored resample matrix, column by column
  expect_equal(rs$value, as.numeric(boot$M_G2))
  lean <- bootstrap_pathways(context_tree(seqs, max_depth = 1L), iter = 5L,
                             keep_resamples = FALSE)
  expect_error(as.data.frame(lean, what = "resamples"),
               class = "transitiontrees_no_data")
})

test_that("plot methods for bootstrap results draw without messages", {
  set.seed(1)
  seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  boot <- bootstrap_pathways(context_tree(seqs, max_depth = 1L), iter = 10L)
  p <- plot(boot)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  lean <- bootstrap_pathways(context_tree(seqs, max_depth = 1L), iter = 5L,
                             keep_resamples = FALSE)
  expect_error(plot_pathway_resamples(lean),
               class = "transitiontrees_no_data")
})
