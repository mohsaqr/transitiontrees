# context_tree(group = ...), model_fit(), n_nodes(), summary() and as.data.frame() of a group, compare_groups():
# one tree per group over a shared alphabet. Each group's tree must be PST::pstree() of that group's sequences.
# compare_groups() then tests every context of any group's tree: the counts of the context in each group
# (PST::cprob(prob = FALSE) of the group's sequences), the Jensen-Shannon divergence of the groups' next-state
# distributions in bits (entropy::entropy.plugin() of the mixture minus the weighted entropies), and the G2 of the
# context's usage (the group x [this context, the other tested contexts of its depth] table; stats::glm() Poisson
# deviance of independence). The permutation p-values are drawn inside the function with no argument to supply
# the permutations, and are not compared.

grouped_input <- function() {
  sequences <- ref_chain(150L, 20L, 3L, 51L)
  # The second group from a different chain, so that the groups differ.
  sequences[seq(2L, 150L, by = 2L), ] <- ref_chain(75L, 20L, 3L, 52L)
  list(sequences = sequences, group = rep(c("odd", "even"), length.out = 150L), alphabet = LETTERS[1:3],
       max_depth = 3L, min_count = 5L, ymin = 0.001)
}

fit_groups <- function(input) {
  context_tree(input$sequences, group = input$group, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

pst_groups <- function(input) {
  lapply(ref_split_rows(input$sequences, input$group),
         \(rows) pst_fit(rows, input$alphabet, input$max_depth, input$min_count, input$ymin))
}

test_that("each group's tree equals PST::pstree() of the group's sequences", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- grouped_input()
  groups <- fit_groups(input)
  psts <- pst_groups(input)
  theirs <- do.call(rbind, lapply(names(psts), \(group) {
    nodes <- pst_node_table(psts[[group]])
    probs <- pst_prob_list(psts[[group]])
    cbind(group = group, nodes, next_probability = vapply(nodes$pathway, \(context) max(probs[[context]]), numeric(1)))
  }))
  expect_reference(subset(as.data.frame(groups), select = c(group, pathway, depth, count, next_probability)), theirs,
                   reference = "PST", covers = c("context_tree()", "as.data.frame.transitiontrees_group()"),
                   quantity = "per group: context depth, count and modal probability", tolerance = 1e-10,
                   key = c("group", "pathway"))
  expect_reference(n_nodes(groups), vapply(psts, pst_n_nodes, numeric(1)), reference = "PST", covers = "n_nodes()",
                   quantity = "nodes per group", tolerance = 0)
})

test_that("model_fit() and summary() of a group report each group's PST logLik, df and nobs", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- grouped_input()
  groups <- fit_groups(input)
  psts <- pst_groups(input)
  theirs <- do.call(rbind, lapply(names(psts), \(group) {
    pst <- psts[[group]]
    data.frame(group = group, logLik = as.numeric(pst_logLik(pst)), df = pst_freepar(pst), nobs = pst_nobs(pst))
  }))
  expect_reference(subset(model_fit(groups), select = c(group, logLik, df, nobs)), theirs, reference = "PST",
                   covers = "model_fit()", quantity = "per group: logLik, df and nobs", tolerance = 1e-10,
                   measure = "relative", key = "group")
  information <- transform(theirs, AIC = -2 * logLik + 2 * df, BIC = -2 * logLik + log(nobs) * df, perplexity = exp(-logLik / nobs))
  expect_reference(subset(model_fit(groups), select = c(group, AIC, BIC, perplexity)),
                   subset(information, select = c(group, AIC, BIC, perplexity)),
                   reference = "independent: AIC, BIC and exp(-logLik / nobs) from each group's PST logLik, df and nobs",
                   covers = "model_fit()", quantity = "per group: AIC, BIC and perplexity", tolerance = 1e-10,
                   measure = "relative", key = "group")
  sizes <- do.call(rbind, lapply(names(psts), \(group) {
    data.frame(group = group, n_seq = sum(input$group == group), n_obs = pst_nobs(psts[[group]]),
               n_nodes = pst_n_nodes(psts[[group]]), logLik = as.numeric(pst_logLik(psts[[group]])))
  }))
  expect_reference(subset(as.data.frame(summary(groups)), select = c(group, n_seq, n_obs, n_nodes, logLik)), sizes,
                   reference = "PST", covers = "summary.transitiontrees_group()",
                   quantity = "per group: sequences, observations, nodes and logLik", tolerance = 1e-10,
                   measure = "relative", key = "group")
})

# The reference side of compare_groups(): every context of any group's PST tree.
ref_group_comparison <- function(input, psts) {
  per_group <- ref_split_rows(input$sequences, input$group)
  contexts <- unique(unlist(lapply(psts, \(pst) pst_node_table(pst)$pathway)))
  depth <- vapply(contexts, \(context) if (context == "(start)") 0L else length(strsplit(context, " -> ", fixed = TRUE)[[1L]]),
                  integer(1), USE.NAMES = FALSE)
  # Count of each context and next state in each group: PST::cprob() of the group's sequences, all contexts.
  tables <- lapply(per_group, \(rows) {
    states <- pst_sequences(rows, input$alphabet)
    lapply(stats::setNames(sort(unique(depth)), sort(unique(depth))), \(d) {
      counts <- suppressMessages(PST::cprob(states, L = d, prob = FALSE))
      rownames(counts) <- pst_to_arrow(if (d == 0L) "e" else rownames(counts))
      counts[, input$alphabet, drop = FALSE]
    })
  })
  count_matrix <- \(context, d) t(vapply(names(per_group), \(group) {
    table <- tables[[group]][[as.character(d)]]
    if (context %in% rownames(table)) as.numeric(table[context, ]) else rep(0, length(input$alphabet))
  }, numeric(length(input$alphabet))))
  matrices <- Map(count_matrix, contexts, depth)
  jsd <- vapply(matrices, \(counts) {
    totals <- rowSums(counts)
    seen <- totals > 0
    weights <- totals[seen] / sum(totals[seen])
    distributions <- counts[seen, , drop = FALSE] / totals[seen]
    mixture <- colSums(weights * distributions)
    each <- vapply(seq_len(nrow(distributions)), \(i) ref_entropy(distributions[i, ], unit = "log2"), numeric(1))
    max(ref_entropy(mixture, unit = "log2") - sum(weights * each), 0)
  }, numeric(1))
  depth_totals <- lapply(stats::setNames(sort(unique(depth)), sort(unique(depth))), \(d)
    Reduce(`+`, lapply(matrices[depth == d], rowSums)))
  usage <- vapply(seq_along(contexts), \(i) {
    used <- rowSums(matrices[[i]])
    others <- pmax(depth_totals[[as.character(depth[[i]])]] - used, 0)
    cells <- data.frame(count = c(used, others), group = rep(names(used), 2L), side = rep(c("this", "other"), each = length(used)))
    # The root is the only context of depth 0: its table has an empty column and a G2 of 0.
    if (depth[[i]] == 0L) 0 else stats::glm(count ~ group + side, family = stats::poisson(), data = cells)$deviance
  }, numeric(1))
  counts <- vapply(matrices, rowSums, numeric(length(per_group)))
  rows <- data.frame(pathway = contexts, depth = depth, count_total = colSums(counts), stringsAsFactors = FALSE)
  rows <- cbind(rows, stats::setNames(as.data.frame(t(counts)), paste0("count_", names(per_group))))
  rows$jsd_bits <- unname(jsd)
  rows$usage_g2 <- ifelse(contexts == "(start)", NA_real_, usage)
  list(table = rows, behavioral = sum(colSums(counts) * jsd), usage = sum(usage))
}

tree_distance_reference <- function(psts) {
  pst_a <- psts[["odd"]]
  pst_b <- psts[["even"]]
  nodes_a <- pst_node_table(pst_a)
  nodes_b <- pst_node_table(pst_b)
  contexts <- unique(c(nodes_a$pathway, nodes_b$pathway))
  count_of <- \(nodes) ifelse(contexts %in% nodes$pathway, nodes$count[match(contexts, nodes$pathway)], 0)
  symmetric <- vapply(contexts, \(context) {
    p <- pst_query(pst_a, context)
    q <- pst_query(pst_b, context)
    (ref_kl(p, q, unit = "log") + ref_kl(q, p, unit = "log")) / 2
  }, numeric(1))
  weight <- count_of(nodes_a) + count_of(nodes_b)
  sum(symmetric * weight) / sum(weight)
}

test_that("compare_groups() gives each context's group counts, JSD and usage G2, and the omnibus statistics", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  input <- grouped_input()
  groups <- fit_groups(input)
  comparison <- compare_groups(groups, iter = 19, seed = 1)
  theirs <- ref_group_comparison(input, pst_groups(input))
  ours <- subset(as.data.frame(comparison), select = c(pathway, depth, count_total, count_odd, count_even, jsd_bits, usage_g2))
  # The usage G2 is a glm deviance, converged by IRLS to its default epsilon 1e-8: tolerance 1e-6.
  expect_reference(subset(ours, select = -usage_g2), subset(theirs$table, select = -usage_g2), reference = "entropy",
                   covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                   quantity = "per context: group counts and Jensen-Shannon divergence (bits)",
                   tolerance = 1e-12, key = "pathway")
  expect_reference(subset(ours, select = c(pathway, usage_g2)), subset(theirs$table, select = c(pathway, usage_g2)),
                   reference = "stats", covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                   quantity = "per context: usage G2",
                   tolerance = 1e-6, key = "pathway")
  overview <- ref_tag(subset(as.data.frame(summary(comparison)), select = c(behavioral_jsd, usage_g2)), test = "omnibus")
  expect_reference(overview, data.frame(test = "omnibus", behavioral_jsd = theirs$behavioral, usage_g2 = theirs$usage),
                   reference = "stats", key = "test",
                   covers = c("compare_groups()", "summary.transitiontrees_group_comparison()"),
                   quantity = "omnibus: count-weighted JSD and summed usage G2", tolerance = 1e-6, measure = "relative")
  distance <- subset(as.data.frame(comparison, what = "distance"), select = c(group_a, group_b, distance))
  expect_reference(distance, data.frame(group_a = "even", group_b = "odd",
                                        distance = tree_distance_reference(pst_groups(input))),
                   reference = "independent: count-weighted mean KL over the union of contexts, of PST::query() distributions (entropy::KL.plugin)",
                   covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                   quantity = "distance between the groups' trees", tolerance = 1e-12, key = c("group_a", "group_b"))
})
