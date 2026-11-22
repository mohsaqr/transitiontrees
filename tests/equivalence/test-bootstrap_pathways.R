# bootstrap_pathways(), summary() and as.data.frame() of its result: the case bootstrap of every pathway of a tree.
# The resampled indices are drawn inside the function (seed = only), so the bootstrap part is checked by seeded Monte
# Carlo, fixed before the run:
#   1. against boot::boot() (Canty & Ripley), the established case-resampling implementation, run with its own seed
#      on a statistic written in base R from the definitions in ?bootstrap_pathways (each pathway's count, modal
#      next-state probability, KL divergence from its parent in bits, G2 against the parent, the band test of the
#      count, the share of resamples whose G2 exceeds qchisq(0.95, k - 1), the share keeping the original modal
#      flip). The two runs are independent Monte Carlo estimates of the same bootstrap expectation. Means are compared
#      in units of the standard error of their difference, sd * sqrt(1 / B1 + 1 / B2) with sd the bootstrap SD from
#      boot's replicates, within 4; rates within the worst-case bound 4 * sqrt(0.25 * (1 / B1 + 1 / B2)).
#   2. against the exact bootstrap moments of a count (a sum over resampled sequences): mean = the observed count,
#      variance = n * the population variance of the per-sequence counts (Efron & Tibshirani 1993, ch. 5), within
#      4 Monte Carlo standard errors of the mean and of the SD (the SD's from the kurtosis of the sum).
# The pathways' observed statistics (no resampling) are compared with PST::cprob() counts of
# the same sequences: counts and modal next states, entropy::KL.plugin() of the distribution and its parent's.
# Sequences: simulab::simulate_markov() chains of a known three-state chain, cut to lengths 10 to 25.

alphabet <- c("A", "B", "C")

boot_sequences <- function(seed) {
  transition <- matrix(c(0.6, 0.3, 0.1,  0.2, 0.5, 0.3,  0.35, 0.15, 0.5), nrow = 3L, byrow = TRUE)
  chains <- simulab::simulate_markov(n = 200L, transition = transition, chain_length = 25L, states = alphabet,
                                     initial = rep(1 / 3, 3L), seed = seed)
  wide <- as.data.frame(chains, what = "wide")
  sequences <- as.matrix(wide[, paste0("S", seq_len(25L))])
  keep <- 10L + (seq_len(nrow(sequences)) %% 16L)
  sequences[col(sequences) > keep] <- NA_character_
  dimnames(sequences) <- list(NULL, paste0("T", seq_len(25L)))
  sequences
}

# Per-sequence next-state counts of every pathway: an array sequence x pathway x next state. A pathway "A -> B" is
# the history A then B (oldest first); "(start)" counts every position.
ref_sequence_counts <- function(sequences, pathways) {
  rows <- lapply(seq_len(nrow(sequences)), \(s) {
    states <- stats::na.omit(sequences[s, ])
    depths <- lengths(strsplit(setdiff(pathways, "(start)"), " -> ", fixed = TRUE))
    windows <- lapply(sort(unique(depths)), \(d) {
      starts <- seq_len(max(length(states) - d, 0L))
      data.frame(context = vapply(starts, \(t) paste(states[t:(t + d - 1L)], collapse = " -> "), character(1)),
                 next_state = states[starts + d], stringsAsFactors = FALSE)
    })
    rbind(data.frame(sequence = s, context = "(start)", next_state = states),
          cbind(sequence = s, do.call(rbind, windows)))
  })
  long <- do.call(rbind, rows)
  table(factor(long$sequence, levels = seq_len(nrow(sequences))), factor(long$context, levels = pathways),
        factor(long$next_state, levels = alphabet))
}

# The statistics of every pathway from its pooled count vector and its parent's, as bootstrap_pathways() defines them.
ref_pathway_stats <- function(pooled, pathways) {
  parent <- ref_parent(pathways)
  k <- length(alphabet)
  one <- function(j) {
    counts <- pooled[j, ]
    total <- sum(counts)
    out <- c(count = total, prob = NA_real_, kl = NA_real_, g2 = NA_real_, flip = NA_real_)
    if (total == 0 || pathways[[j]] == "(start)") {
      if (total > 0) out[["prob"]] <- max(counts) / total
      return(out)
    }
    parent_counts <- pooled[match(parent[[j]], pathways), ]
    p <- counts / total
    q <- parent_counts / sum(parent_counts)
    seen <- counts > 0
    out[["prob"]] <- max(counts) / total
    out[["kl"]] <- sum(p[seen] * log2(p[seen] / q[seen]))
    out[["g2"]] <- 2 * sum(counts[seen] * log(counts[seen] / (total * q[seen])))
    out[["flip"]] <- as.numeric(which.max(counts) != which.max(parent_counts))
    out
  }
  stats <- vapply(seq_along(pathways), one, numeric(5))
  stats
}

# boot::boot() of the resampled statistics: per resample, per pathway, the five statistics and the band and G2 flags.
ref_boot <- function(per_sequence, pathways, original, replicates, seed) {
  critical <- stats::qchisq(0.95, df = length(alphabet) - 1L)
  statistic <- function(data, index) {
    pooled <- apply(per_sequence[index, , , drop = FALSE], c(2L, 3L), sum)
    s <- ref_pathway_stats(pooled, pathways)
    band <- as.numeric(s["count", ] >= 0.5 * original["count", ] & s["count", ] <= 1.5 * original["count", ])
    informative <- as.numeric(!is.na(s["g2", ]) & s["g2", ] > critical)
    keeps_flip <- as.numeric(s["flip", ] == original["flip", ])
    c(s["count", ], s["prob", ], s["kl", ], s["g2", ], band, informative, keeps_flip)
  }
  set.seed(seed)
  result <- boot::boot(data = seq_len(dim(per_sequence)[[1L]]), statistic = statistic, R = replicates)
  block <- \(i) result$t[, (i - 1L) * length(pathways) + seq_along(pathways), drop = FALSE]
  column_mean <- \(m) apply(m, 2L, \(v) mean(v, na.rm = TRUE))
  column_sd <- \(m) apply(m, 2L, \(v) stats::sd(v, na.rm = TRUE))
  names(pathways) <- NULL
  data.frame(pathway = pathways,
             mean_count = column_mean(block(1L)), sd_count = column_sd(block(1L)),
             mean_next_probability = column_mean(block(2L)), sd_next_probability = column_sd(block(2L)),
             mean_divergence = column_mean(block(3L)), sd_divergence = column_sd(block(3L)),
             mean_G2 = column_mean(block(4L)), sd_G2 = column_sd(block(4L)),
             stability_rate = column_mean(block(5L)),
             informative_rate = ifelse(pathways == "(start)", NA_real_, colSums(block(6L)) / replicates),
             flip_consistency = column_mean(block(7L)), stringsAsFactors = FALSE)
}

iterations <- 2000L

boot_case <- function(seed) {
  sequences <- boot_sequences(seed)
  tree <- context_tree(sequences, max_depth = 3L, min_count = 5L, alphabet = alphabet)
  found <- summary(bootstrap_pathways(tree, iter = iterations, seed = seed))
  pathways <- found$pathway
  per_sequence <- ref_sequence_counts(sequences, pathways)
  original <- ref_pathway_stats(apply(per_sequence, c(2L, 3L), sum), pathways)
  list(sequences = sequences, found = found, pathways = pathways, per_sequence = per_sequence, original = original,
       boot = ref_boot(per_sequence, pathways, original, iterations, seed + 1000L))
}

cases <- lapply(stats::setNames(c(1L, 2L, 3L), paste("seed", c(1L, 2L, 3L))), boot_case)

test_that("the bootstrap means and rates of every pathway equal boot::boot() of the definition (seeded Monte Carlo)", {
  skip_without_reference("boot")
  invisible(lapply(names(cases), \(name) {
    case <- cases[[name]]
    reference <- case$boot
    found <- case$found
    found_rows <- match(reference$pathway, found$pathway)
    se_factor <- sqrt(2 / iterations)
    means <- c("count", "next_probability", "divergence", "G2")
    standardise <- \(table, rows) {
      out <- data.frame(pathway = reference$pathway)
      out[paste0("mean_", means)] <- lapply(means, \(m) table[[paste0("mean_", m)]][rows] /
                                              (reference[[paste0("sd_", m)]] * se_factor))
      out
    }
    expect_reference(standardise(found, found_rows), standardise(reference, seq_len(nrow(reference))),
                     reference = "boot", covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: bootstrap means in standard errors of the difference", name), tolerance = 4,
                     key = "pathway",
                     note = "both sides divided by sd * sqrt(1/2000 + 1/2000), sd from boot's replicates; 4 SE fixed before the run")
    # Rates are shares of 2000 resamples: worst-case SE of a difference sqrt(0.25 * 2 / 2000), 4 of them.
    expect_reference(found, subset(reference, select = c(pathway, stability_rate, informative_rate, flip_consistency)),
                     reference = "boot", covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: stability, informative and flip-consistency rates", name),
                     tolerance = 4 * sqrt(0.25 * 2 / iterations), key = "pathway")
  }))
})

test_that("as.data.frame() of the bootstrap is its summary, the same boot::boot() means", {
  skip_without_reference("boot")
  case <- cases[[1L]]
  found <- as.data.frame(bootstrap_pathways(context_tree(case$sequences, max_depth = 3L, min_count = 5L,
                                                         alphabet = alphabet), iter = iterations, seed = 1L))
  expect_reference(found, subset(case$boot, select = c(pathway, stability_rate, informative_rate, flip_consistency)),
                   reference = "boot", covers = c("bootstrap_pathways()", "as.data.frame.transitiontrees_bootstrap()"),
                   quantity = "seed 1: rates from as.data.frame()", tolerance = 4 * sqrt(0.25 * 2 / iterations),
                   key = "pathway")
})

test_that("the bootstrap mean and SD of every count equal the exact bootstrap moments of a sum (seeded Monte Carlo)", {
  invisible(lapply(names(cases), \(name) {
    case <- cases[[name]]
    per_count <- apply(case$per_sequence, c(1L, 2L), sum)
    n <- nrow(per_count)
    centred <- sweep(per_count, 2L, colMeans(per_count))
    m2 <- colMeans(centred^2)
    m4 <- colMeans(centred^4)
    sigma <- sqrt(n * m2)
    kurtosis_of_sum <- 3 + (m4 / m2^2 - 3) / n
    rows <- match(case$pathways, case$found$pathway)
    se_mean <- sigma / sqrt(iterations)
    se_sd <- sigma * sqrt((kurtosis_of_sum - 1) / (4 * iterations))
    found <- data.frame(pathway = case$pathways, mean_count = case$found$mean_count[rows] / se_mean,
                        sd_count = case$found$sd_count[rows] / se_sd)
    expected <- data.frame(pathway = case$pathways, mean_count = colSums(per_count) / se_mean, sd_count = sigma / se_sd)
    expect_reference(found, expected,
                     reference = "independent: exact bootstrap mean and variance of a sum of resampled sequence counts",
                     covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: mean and SD of the count in Monte Carlo standard errors", name),
                     tolerance = 4, key = "pathway",
                     note = "SE of the mean sigma / sqrt(B); of the SD sigma * sqrt((kurtosis - 1) / (4 B)); 4 SE fixed before the run")
  }))
})

test_that("the observed statistics of every pathway equal PST::cprob() counts and entropy::KL.plugin()", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  invisible(lapply(names(cases), \(name) {
    case <- cases[[name]]
    # PST::cprob(prob = FALSE): the next-state counts of every context of each length (pstree() itself stops on
    # these ragged sequences, in its likelihood step, TraMineR::seqdecomp()).
    states <- pst_sequences(case$sequences, alphabet)
    tables <- lapply(0:3, \(d) {
      counts <- suppressMessages(PST::cprob(states, L = d, prob = FALSE))
      rownames(counts) <- pst_to_arrow(if (d == 0L) "e" else rownames(counts))
      counts[, alphabet, drop = FALSE]
    })
    all_counts <- do.call(rbind, tables)
    distribution <- \(pathway) {
      counts <- as.numeric(all_counts[pathway, ])
      stats::setNames(counts / sum(counts), alphabet)
    }
    expected <- do.call(rbind, lapply(case$pathways, \(pathway) {
      p <- distribution(pathway)
      data.frame(pathway = pathway, count = sum(all_counts[pathway, ]), likely_next = names(p)[which.max(p)],
                 next_probability = max(p),
                 divergence = if (pathway == "(start)") NA_real_ else ref_kl(p, distribution(ref_parent(pathway))),
                 stringsAsFactors = FALSE)
    }))
    expect_reference(case$found, subset(expected, select = c(pathway, count, likely_next, next_probability)),
                     reference = "PST", covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: observed count and modal next state", name), tolerance = 1e-12,
                     key = "pathway", note = "PST::cprob(prob = FALSE) of the same sequences")
    expect_reference(case$found, subset(expected, select = c(pathway, divergence)), reference = "entropy",
                     covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: observed divergence from the parent (bits)", name), tolerance = 1e-12,
                     key = "pathway", note = "entropy::KL.plugin() of PST::cprob() distributions of the pathway and its parent")
  }))
})
