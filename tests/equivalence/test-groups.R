# context_tree(group = ...), model_fit(), n_nodes(), summary() and as.data.frame() of a group, compare_groups():
# one tree per group over a shared alphabet. Each group's tree must be PST::pstree() of that group's sequences.
# compare_groups() then tests every context of any group's tree: the counts of the context in each group
# (PST::cprob(prob = FALSE) of the group's sequences), the Jensen-Shannon divergence of the groups' next-state
# distributions in bits (entropy::entropy.plugin() of the mixture minus the weighted entropies), and the G2 of the
# context's usage (the group x [this context, the other tested contexts of its depth] table; stats::glm() Poisson
# deviance of independence). The permutation p-values are drawn inside the function with no argument to supply
# the permutations: they are compared by seeded Monte Carlo with a label permutation written out in base R from
# ?compare_groups (the same contexts, statistics and tie rule, its own seed), cell by cell in standard errors of
# the difference of two Monte Carlo p-values, sqrt(p (1 - p) (1 / B1 + 1 / B2)) with p the mean of the two, within
# 4 (fixed before the run), over three seeds.

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
  # The behavioral omnibus is a sum of entropy::entropy.plugin() divergences; the usage omnibus a sum of glm deviances.
  expect_reference(subset(overview, select = c(test, behavioral_jsd)),
                   data.frame(test = "omnibus", behavioral_jsd = theirs$behavioral),
                   reference = "entropy", key = "test",
                   covers = c("compare_groups()", "summary.transitiontrees_group_comparison()"),
                   quantity = "omnibus: count-weighted JSD", tolerance = 1e-12, measure = "relative")
  expect_reference(subset(overview, select = c(test, usage_g2)), data.frame(test = "omnibus", usage_g2 = theirs$usage),
                   reference = "stats", key = "test",
                   covers = c("compare_groups()", "summary.transitiontrees_group_comparison()"),
                   quantity = "omnibus: summed usage G2", tolerance = 1e-6, measure = "relative")
  distance <- subset(as.data.frame(comparison, what = "distance"), select = c(group_a, group_b, distance))
  expect_reference(distance, data.frame(group_a = "even", group_b = "odd",
                                        distance = tree_distance_reference(pst_groups(input))),
                   reference = "independent: count-weighted mean KL over the union of contexts, of PST::query() distributions (entropy::KL.plugin)",
                   covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                   quantity = "distance between the groups' trees", tolerance = 1e-12, key = c("group_a", "group_b"))
})

# The permutation null of compare_groups() written out: sequences keep their counts, labels are shuffled
# (base::sample() of the label vector, the reference's own seed), and every context's JSD and usage G2 and the two
# omnibus totals are recomputed; p = (1 + #{permuted >= observed - 1e-12}) / (B + 1), as ?compare_groups defines.
ref_group_permutation <- function(input, contexts, iterations, seed) {
  per_sequence <- ref_sequence_counts(input$sequences, contexts, input$alphabet)
  dims <- dim(per_sequence)
  flat <- matrix(per_sequence, nrow = dims[[1L]])
  depth <- vapply(contexts, \(context) if (context == "(start)") 0L else length(strsplit(context, " -> ", fixed = TRUE)[[1L]]),
                  integer(1), USE.NAMES = FALSE)
  entropy_bits <- \(p) -sum(p[p > 0] * log2(p[p > 0]))
  jsd_of <- \(counts) {
    totals <- rowSums(counts)
    seen <- totals > 0
    if (sum(seen) < 2L) return(0)
    weights <- totals[seen] / sum(totals[seen])
    distributions <- counts[seen, , drop = FALSE] / totals[seen]
    max(entropy_bits(colSums(weights * distributions)) -
          sum(weights * apply(distributions, 1L, entropy_bits)), 0)
  }
  g2_of <- \(table) {
    expected <- outer(rowSums(table), colSums(table)) / sum(table)
    seen <- table > 0 & expected > 0
    if (sum(table) == 0) 0 else 2 * sum(table[seen] * log(table[seen] / expected[seen]))
  }
  statistics <- function(labels) {
    indicator <- stats::model.matrix(~ 0 + factor(labels))
    pooled <- array(crossprod(indicator, flat), c(ncol(indicator), dims[[2L]], dims[[3L]]))
    used <- apply(pooled, c(1L, 2L), sum)
    totals <- vapply(seq_along(contexts), \(i) rowSums(used[, depth == depth[[i]], drop = FALSE]), numeric(nrow(used)))
    jsd <- vapply(seq_along(contexts), \(i) jsd_of(pooled[, i, ]), numeric(1))
    usage <- vapply(seq_along(contexts), \(i) g2_of(cbind(used[, i], pmax(totals[, i] - used[, i], 0))), numeric(1))
    list(jsd = jsd, usage = usage, omnibus = c(behavioral = sum(colSums(used) * jsd), usage = sum(usage)))
  }
  observed <- statistics(input$group)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  permuted <- replicate(iterations, statistics(sample(input$group)), simplify = FALSE)
  exceed <- \(part) rowSums(vapply(permuted, \(s) s[[part]] >= observed[[part]] - 1e-12, logical(length(observed[[part]]))))
  p_of <- \(part) (1 + exceed(part)) / (iterations + 1)
  list(contexts = data.frame(pathway = contexts, jsd_p = p_of("jsd"),
                             usage_p = ifelse(contexts == "(start)", NA_real_, p_of("usage")), stringsAsFactors = FALSE),
       omnibus = data.frame(axis = c("behavioral", "usage"), p_value = unname(p_of("omnibus")), stringsAsFactors = FALSE))
}

test_that("compare_groups() permutation p-values agree with a label permutation written out (seeded Monte Carlo)", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- grouped_input()
  groups <- fit_groups(input)
  contexts <- unique(unlist(lapply(pst_groups(input), \(pst) pst_node_table(pst)$pathway)))
  iterations <- 999L
  scaled <- function(found, theirs, key, columns) {
    rows <- match(theirs[[key]], found[[key]])
    se <- lapply(stats::setNames(columns, columns), \(column) {
      # The mean of the two, kept off 0 and 1 so two p-values of exactly 1 stay equal, not 0 / 0.
      p <- pmin(pmax((found[[column]][rows] + theirs[[column]]) / 2, 0.5 / iterations), 1 - 0.5 / iterations)
      sqrt(p * (1 - p) * 2 / iterations)
    })
    side <- \(table, rows) cbind(stats::setNames(data.frame(theirs[[key]], stringsAsFactors = FALSE), key),
                                 as.data.frame(lapply(stats::setNames(columns, columns), \(column) table[[column]][rows] / se[[column]])))
    list(found = side(found, rows), theirs = side(theirs, seq_len(nrow(theirs))))
  }
  invisible(lapply(c(1L, 2L, 3L), \(seed) {
    comparison <- compare_groups(groups, iter = iterations, seed = seed)
    theirs <- ref_group_permutation(input, contexts, iterations, seed + 1000L)
    per_context <- scaled(as.data.frame(comparison), theirs$contexts, "pathway", c("jsd_p", "usage_p"))
    expect_reference(per_context$found, per_context$theirs,
                     reference = "independent: label permutation of ?compare_groups written out in base R",
                     covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                     quantity = sprintf("seed %d: per-context permutation p-values in standard errors", seed),
                     tolerance = 4, key = "pathway",
                     note = "999 permutations each side; both divided by sqrt(p (1 - p) (1/999 + 1/999)), p the mean of the two; 4 SE fixed before the run")
    omnibus <- scaled(as.data.frame(comparison, what = "omnibus"), theirs$omnibus, "axis", "p_value")
    expect_reference(omnibus$found, omnibus$theirs,
                     reference = "independent: label permutation of ?compare_groups written out in base R",
                     covers = c("compare_groups()", "as.data.frame.transitiontrees_group_comparison()"),
                     quantity = sprintf("seed %d: omnibus permutation p-values in standard errors", seed),
                     tolerance = 4, key = "axis",
                     note = "999 permutations each side; both divided by sqrt(p (1 - p) (1/999 + 1/999)), p the mean of the two; 4 SE fixed before the run")
  }))
})
