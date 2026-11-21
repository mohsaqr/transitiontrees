# tree_distance(), compare_trees(): how far two trees' predictions are apart. Over the union of the two trees'
# contexts, each tree's next-state distribution (the context's own, or that of its longest suffix in the tree:
# PST::query()), the Kullback-Leibler divergence each way in nats (entropy::KL.plugin()), their mean, weighted by
# the two counts of the context (0 where a tree lacks it). PST::pdist() is a different quantity (a Monte Carlo
# log-likelihood ratio of generated sequences) and is not used. The permutation null of compare_trees() is drawn
# inside the function with no argument to supply the permutations: only its observed part is compared here.

tree_pair <- function(a, b, alphabet, max_depth = 3L, min_count = 5L, ymin = 0.001) {
  list(a = context_tree(a, max_depth = max_depth, min_count = min_count, smoothing = list("floor", ymin = ymin), alphabet = alphabet),
       b = context_tree(b, max_depth = max_depth, min_count = min_count, smoothing = list("floor", ymin = ymin), alphabet = alphabet),
       pst_a = pst_fit(a, alphabet, max_depth, min_count, ymin), pst_b = pst_fit(b, alphabet, max_depth, min_count, ymin))
}

pairs <- function() {
  chain <- ref_chain(160L, 20L, 3L, 41L)
  list(halves = tree_pair(chain[1:80, ], chain[81:160, ], LETTERS[1:3]),
       two_chains = tree_pair(ref_chain(100L, 20L, 3L, 42L), ref_chain(100L, 20L, 3L, 43L), LETTERS[1:3]),
       four_states = tree_pair(ref_chain(120L, 25L, 4L, 44L), ref_chain(90L, 25L, 4L, 45L), LETTERS[1:4], max_depth = 2L))
}

# The breakdown per context, written out on PST's numbers.
ref_breakdown <- function(pst_a, pst_b) {
  nodes_a <- pst_node_table(pst_a)
  nodes_b <- pst_node_table(pst_b)
  contexts <- unique(c(nodes_a$pathway, nodes_b$pathway))
  count_of <- \(nodes) ifelse(contexts %in% nodes$pathway, nodes$count[match(contexts, nodes$pathway)], 0)
  kl <- \(from, to) vapply(contexts, \(context) ref_kl(pst_query(from, context), pst_query(to, context), unit = "log"), numeric(1))
  ab <- kl(pst_a, pst_b)
  ba <- kl(pst_b, pst_a)
  data.frame(pathway = contexts, count_a = count_of(nodes_a), count_b = count_of(nodes_b), divergence_ab = unname(ab),
             divergence_ba = unname(ba), divergence_sym = unname((ab + ba) / 2), stringsAsFactors = FALSE)
}

ref_distance <- function(breakdown, symmetric = TRUE) {
  per <- if (symmetric) breakdown$divergence_sym else breakdown$divergence_ab
  weight <- breakdown$count_a + breakdown$count_b
  sum(per * weight) / sum(weight)
}

reference <- "independent: count-weighted mean KL over the union of contexts, of PST::query() distributions (entropy::KL.plugin)"

test_that("tree_distance() is the count-weighted divergence of PST's distributions over the union of contexts", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  pairs <- pairs()
  ours <- vapply(pairs, \(pair) c(symmetric = tree_distance(pair$a, pair$b), directed = tree_distance(pair$a, pair$b, symmetric = FALSE)),
                 numeric(2))
  theirs <- vapply(pairs, \(pair) {
    breakdown <- ref_breakdown(pair$pst_a, pair$pst_b)
    c(symmetric = ref_distance(breakdown), directed = ref_distance(breakdown, symmetric = FALSE))
  }, numeric(2))
  expect_reference(ours, theirs, reference = reference, covers = "tree_distance()", quantity = "distance between two trees",
                   tolerance = 1e-12)
})

test_that("compare_trees() gives the observed distance and the per-context breakdown", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  pairs <- pairs()
  results <- lapply(pairs, \(pair) compare_trees(pair$a, pair$b, iter = 19, seed = 1))
  ours <- do.call(rbind, lapply(names(results), \(name) cbind(pair = name, as.data.frame(results[[name]]))))
  theirs <- do.call(rbind, lapply(names(pairs), \(name) cbind(pair = name, ref_breakdown(pairs[[name]]$pst_a, pairs[[name]]$pst_b))))
  expect_reference(ours, theirs, reference = reference, covers = c("compare_trees()", "as.data.frame.transitiontrees_comparison()"),
                   quantity = "per-context counts and divergences", tolerance = 1e-12, key = c("pair", "pathway"))
  overview <- do.call(rbind, lapply(names(results), \(name) cbind(pair = name, subset(summary(results[[name]]), select = c(distance, n_pathways)))))
  expected <- do.call(rbind, lapply(names(pairs), \(name) {
    breakdown <- ref_breakdown(pairs[[name]]$pst_a, pairs[[name]]$pst_b)
    data.frame(pair = name, distance = ref_distance(breakdown), n_pathways = nrow(breakdown))
  }))
  expect_reference(overview, expected, reference = reference, covers = c("compare_trees()", "summary.transitiontrees_comparison()"),
                   quantity = "observed distance and number of contexts", tolerance = 1e-12, key = "pair")
})
