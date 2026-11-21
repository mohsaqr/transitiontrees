# prune_tree(), compare_pruning(): removing leaves that do not predict differently from their parent.
# With criterion = "G2", prune_tree() keeps a leaf when G2 = 2 * sum(n_s * log(n_s / (n * q_s))) exceeds
# qchisq(1 - alpha, k - 1), the leaf's counts tested against its parent's distribution q, leaves first, until
# nothing changes. PST::prune(gain = "G2", C) keeps a node when N * KL(p || q) >= C, with p the node's (floor
# smoothed) distribution; C = qchisq(1 - alpha, k - 1) / 2 is the same test. Two conventions differ, both
# declared below: PST's worked example uses C = qchisq(0.95, 1) / 2 (df 1 for every alphabet), and PST tests the
# smoothed child distribution where prune_tree() tests the observed counts, which differ at a context after which
# a state was never observed. (Replaces the G2 checks of local_testing_and_equivalence/test-equiv-PST.R, which
# required a kept-set subset and 85% agreement.)

fit <- function(sequences, alphabet, max_depth, min_count, ymin) {
  context_tree(sequences, max_depth = max_depth, min_count = min_count, smoothing = list("floor", ymin = ymin),
               alphabet = alphabet)
}

# The kept/pruned status of every context of the unpruned tree.
kept_table <- function(all_contexts, kept_contexts) {
  data.frame(pathway = all_contexts, kept = all_contexts %in% kept_contexts, stringsAsFactors = FALSE)
}

ours_kept <- function(tree, ...) {
  kept_table(tree_pathways(tree)$pathway, tree_pathways(prune_tree(tree, ...))$pathway)
}

theirs_kept <- function(pst, C) {
  kept_table(pst_node_table(pst)$pathway, pst_node_table(pst_prune(pst, gain = "G2", C = C))$pathway)
}

# The inputs on which the two agree: the eight configurations with ymin > 0 (PST cannot G2-prune a tree that
# has a zero probability), long chains of 2 to 4 states, the package's engagement and trajectories data.
agreeing_inputs <- function() {
  configs <- Filter(\(input) input$ymin > 0, ref_inputs())
  long_chains <- lapply(list(c(2L, 300L, 3L, 9L), c(3L, 300L, 2L, 19L), c(4L, 400L, 2L, 29L)), \(cfg)
    list(sequences = ref_chain(cfg[[2]], 30L, cfg[[1]], cfg[[4]]), alphabet = LETTERS[seq_len(cfg[[1]])],
         max_depth = cfg[[3]], min_count = 5L, ymin = 0.001))
  names(long_chains) <- c("long2", "long3", "long4")
  engagement_states <- sort(unique(stats::na.omit(unlist(engagement, use.names = FALSE))))
  trajectory_states <- sort(unique(stats::na.omit(as.vector(trajectories))))
  data_sets <- list(
    engagement = list(sequences = as.matrix(engagement), alphabet = engagement_states, max_depth = 4L, min_count = 5L, ymin = 0.001),
    trajectories = list(sequences = trajectories, alphabet = trajectory_states, max_depth = 3L, min_count = 3L, ymin = 0.001))
  c(configs, long_chains, data_sets)
}

# Chains of 4 states (120 sequences of 25, depth 4, min_count 5) on which a leaf after which a state was never
# observed sits between the two tests (seeds found by running both sides on seeds 1 to 9).
sparse_inputs <- function() {
  inputs <- lapply(c(2L, 4L), \(seed) list(sequences = ref_chain(120L, 25L, 4L, seed), alphabet = LETTERS[1:4],
                                           max_depth = 4L, min_count = 5L, ymin = 0.001))
  stats::setNames(inputs, c("seed2", "seed4"))
}

both_sides <- function(inputs, C_of_k) {
  sides <- lapply(names(inputs), \(name) {
    input <- inputs[[name]]
    tree <- fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    pst <- pst_fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    list(ours = cbind(input = name, ours_kept(tree, criterion = "G2", alpha = 0.05)),
         theirs = cbind(input = name, theirs_kept(pst, C = C_of_k(length(input$alphabet)))),
         independent = cbind(input = name, kept_table(pst_node_table(pst)$pathway,
                                                      ref_prune(pst, "G2", cut = stats::qchisq(0.95, length(input$alphabet) - 1L)))))
  })
  lapply(c(ours = "ours", theirs = "theirs", independent = "independent"), \(part) do.call(rbind, lapply(sides, `[[`, part)))
}

matched_C <- \(k) stats::qchisq(0.95, k - 1L) / 2
worked_example_C <- \(k) stats::qchisq(0.95, 1L) / 2

test_that("prune_tree(criterion = 'G2', alpha = 0.05) keeps the nodes PST::prune(gain = 'G2', C = qchisq(0.95, k - 1) / 2) keeps", {
  skip_without_reference(c("PST", "TraMineR"))
  sides <- both_sides(agreeing_inputs(), matched_C)
  expect_reference(sides$ours, sides$theirs, reference = "PST", covers = "prune_tree()",
                   quantity = "context kept by G2 pruning (alpha = 0.05)", tolerance = 0, key = c("input", "pathway"))
  expect_reference(sides$ours, sides$independent, reference = "independent: G2 of each leaf's PST counts against its parent's PST probabilities, leaves first",
                   covers = "prune_tree()", quantity = "context kept by G2 pruning (alpha = 0.05)", tolerance = 0,
                   key = c("input", "pathway"))
})

test_that("PST's worked-example cutoff (df 1 for every alphabet) keeps more contexts", {
  skip_without_reference(c("PST", "TraMineR"))
  sides <- both_sides(agreeing_inputs(), worked_example_C)
  expect_known_difference(sides$ours, sides$theirs, reference = "PST", covers = "prune_tree()",
                          quantity = "context kept: PST cutoff qchisq(0.95, 1) / 2", tolerance = 0,
                          key = c("input", "pathway"),
                          reason = "PST's worked example sets C = qchisq(0.95, 1) / 2, a chi-square on 1 df for any alphabet; prune_tree() tests on k - 1 df, the df of a k-category likelihood-ratio test. With k = 2 the two are the same.",
                          agrees_with = "PST::prune(gain = 'G2', C = qchisq(1 - alpha, k - 1) / 2)")
})

test_that("at a leaf with an unobserved next state, prune_tree() tests the observed counts and PST the smoothed distribution", {
  skip_without_reference(c("PST", "TraMineR"))
  sides <- both_sides(sparse_inputs(), matched_C)
  expect_known_difference(sides$ours, sides$theirs, reference = "PST", covers = "prune_tree()",
                          quantity = "context kept (inputs with a borderline unobserved state)", tolerance = 0,
                          key = c("input", "pathway"),
                          reason = "PST's G2 gain compares the node's floor-smoothed distribution with its parent's; prune_tree() compares the observed counts. They differ only after a context where a state was never seen (smoothing moves it to ymin), by O(ymin), which moves a leaf across the cutoff: seed 4, leaf C -> C -> C -> A; seed 2, the chain C -> D, A -> C -> D, A -> A -> C -> D.",
                          agrees_with = "G2 on the observed counts (the independent computation in the same test)")
  expect_reference(sides$ours, sides$independent, reference = "independent: G2 of each leaf's PST counts against its parent's PST probabilities, leaves first",
                   covers = "prune_tree()", quantity = "context kept (inputs with a borderline unobserved state)",
                   tolerance = 0, key = c("input", "pathway"))
})

test_that("prune_tree() with the KL, AIC and BIC criteria keeps the contexts the criteria keep", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- agreeing_inputs()
  sides <- lapply(names(inputs), \(name) {
    input <- inputs[[name]]
    tree <- fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    pst <- pst_fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    lapply(c("KL", "AIC", "BIC"), \(criterion) list(
      ours = cbind(input = name, criterion = criterion, ours_kept(tree, criterion = criterion, threshold = 0.005)),
      theirs = cbind(input = name, criterion = criterion,
                     kept_table(pst_node_table(pst)$pathway, ref_prune(pst, criterion, cut = 0.005)))))
  })
  flat <- unlist(sides, recursive = FALSE)
  expect_reference(do.call(rbind, lapply(flat, `[[`, "ours")), do.call(rbind, lapply(flat, `[[`, "theirs")),
                   reference = "independent: KL (> 0.005 nats), AIC and BIC leaf tests on PST's counts and probabilities, leaves first",
                   covers = "prune_tree()", quantity = "context kept by the KL, AIC or BIC criterion", tolerance = 0,
                   key = c("input", "criterion", "pathway"))
})

test_that("compare_pruning() counts the nodes each criterion keeps", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- agreeing_inputs()
  ours <- do.call(rbind, lapply(names(inputs), \(name) {
    input <- inputs[[name]]
    tree <- fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    cbind(input = name, compare_pruning(tree))
  }))
  g2 <- do.call(rbind, lapply(names(inputs), \(name) {
    input <- inputs[[name]]
    pst <- pst_fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    data.frame(input = name, criterion = "G2",
               n_nodes = nrow(pst_node_table(pst_prune(pst, gain = "G2", C = matched_C(length(input$alphabet))))))
  }))
  expect_reference(subset(ours, criterion == "G2", select = c(input, criterion, n_nodes)), g2, reference = "PST",
                   covers = "compare_pruning()", quantity = "nodes kept by G2 pruning", tolerance = 0,
                   key = c("input", "criterion"))
  independent <- do.call(rbind, lapply(names(inputs), \(name) {
    input <- inputs[[name]]
    pst <- pst_fit(input$sequences, input$alphabet, input$max_depth, input$min_count, input$ymin)
    n0 <- nrow(pst_node_table(pst))
    kept <- vapply(c("G2", "KL", "AIC", "BIC"), \(criterion) {
      cut <- if (criterion == "G2") stats::qchisq(0.95, length(input$alphabet) - 1L) else 0.005
      length(ref_prune(pst, criterion, cut = cut))
    }, integer(1))
    data.frame(input = name, criterion = names(kept), n_nodes = unname(kept),
               reduction_pct = round(100 * (1 - unname(kept) / n0), 1))
  }))
  expect_reference(ours, independent,
                   reference = "independent: nodes kept by the four leaf tests on PST's counts, and 100 * (1 - kept / all) to one decimal",
                   covers = "compare_pruning()", quantity = "nodes kept and reduction per criterion", tolerance = 1e-12,
                   key = c("input", "criterion"))
})
