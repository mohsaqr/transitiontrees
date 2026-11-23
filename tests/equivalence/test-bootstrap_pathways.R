# bootstrap_pathways(), summary() and as.data.frame() of its result: the case bootstrap of every pathway of a tree.
# The resampled indices are drawn inside the function (seed = only), so the bootstrap part is checked by seeded Monte
# Carlo, fixed before the run:
#   1. against the definitions in ?bootstrap_pathways written in base R (each pathway's count, modal next-state
#      probability, KL divergence from its parent in bits, G2 against the parent, the band test of the count, the
#      share of resamples whose G2 exceeds qchisq(0.95, k - 1), the share keeping the original modal flip),
#      evaluated on case resamples drawn by boot::boot() (Canty & Ripley) with its own seed. boot only draws the
#      indices; the statistics are the test's, so the reference is an independent computation. The two runs are
#      independent Monte Carlo estimates of the same bootstrap distribution, compared cell by cell in units of the
#      Monte Carlo standard error of their difference, within 4 (fixed before the run):
#      - means: sd * sqrt(1 / B1 + 1 / B2), sd the bootstrap SD from the reference's replicates;
#      - SDs: sd * sqrt((kurtosis - 1) / 4 * (1 / B1 + 1 / B2)), kurtosis of the reference's replicates;
#      - rates and p_stability: sqrt(r (1 - r) (1 / B1 + 1 / B2)), r the mean of the two estimates (kept within
#        [0.5 / B, 1 - 0.5 / B] so that two rates of exactly 0 or 1 compare as equal, not as 0 / 0);
#      - percentile CI bounds: on the probability scale, the reference's empirical CDF at the package's bound
#        (the point of [F(x-), F(x)] nearest the nominal level, for ties) against the nominal level 0.025 or 0.975,
#        SE sqrt(a (1 - a) (1 / B1 + 1 / B2)).
#   2. against the exact bootstrap moments of a count (a sum over resampled sequences): mean = the observed count,
#      variance = n * the population variance of the per-sequence counts (Efron & Tibshirani 1993, ch. 5), within
#      4 Monte Carlo standard errors of the mean and of the SD (the SD's from the kurtosis of the sum).
# The pathways' observed statistics (no resampling) are compared with PST::cprob() counts of
# the same sequences: counts and modal next states, entropy::KL.plugin() of the distribution and its parent's, and
# the G2 of the counts against the parent's distribution written out.
# Not compared: the flags `stable` and `informative` (thresholds of the rates, whose Monte Carlo comparison is above).
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
  replicates_of <- list(count = block(1L), next_probability = block(2L), divergence = block(3L), G2 = block(4L))
  out <- data.frame(pathway = pathways,
             mean_count = column_mean(block(1L)), sd_count = column_sd(block(1L)),
             mean_next_probability = column_mean(block(2L)), sd_next_probability = column_sd(block(2L)),
             mean_divergence = column_mean(block(3L)), sd_divergence = column_sd(block(3L)),
             mean_G2 = column_mean(block(4L)), sd_G2 = column_sd(block(4L)),
             stability_rate = column_mean(block(5L)),
             informative_rate = ifelse(pathways == "(start)", NA_real_, colSums(block(6L)) / replicates),
             flip_consistency = column_mean(block(7L)),
             p_stability = (replicates - colSums(block(5L)) + 1) / (replicates + 1), stringsAsFactors = FALSE)
  attr(out, "replicates") <- replicates_of
  out
}

iterations <- 2000L

boot_case <- function(seed) {
  sequences <- boot_sequences(seed)
  tree <- context_tree(sequences, max_depth = 3L, min_count = 5L, alphabet = alphabet)
  found <- summary(bootstrap_pathways(tree, iter = iterations, seed = seed))
  pathways <- found$pathway
  per_sequence <- ref_sequence_counts(sequences, pathways, alphabet)
  original <- ref_pathway_stats(apply(per_sequence, c(2L, 3L), sum), pathways)
  list(sequences = sequences, found = found, pathways = pathways, per_sequence = per_sequence, original = original,
       boot = ref_boot(per_sequence, pathways, original, iterations, seed + 1000L))
}

cases <- lapply(stats::setNames(c(1L, 2L, 3L), paste("seed", c(1L, 2L, 3L))), boot_case)

boot_reference <- "independent: ?bootstrap_pathways definitions on boot::boot() case resamples"
statistics <- c("count", "next_probability", "divergence", "G2")
rate_columns <- c("stability_rate", "informative_rate", "flip_consistency", "p_stability")

# Both sides of a cell divided by the same Monte Carlo standard error of their difference.
standardised_pair <- function(found, reference, columns, se_of) {
  rows <- match(reference$pathway, found$pathway)
  scaled <- lapply(stats::setNames(columns, columns), \(column) {
    se <- se_of(column, found[[column]][rows], reference[[column]])
    list(found = found[[column]][rows] / se, reference = reference[[column]] / se)
  })
  side <- \(part) cbind(data.frame(pathway = reference$pathway, stringsAsFactors = FALSE),
                        as.data.frame(lapply(scaled, `[[`, part)))
  list(found = side("found"), reference = side("reference"))
}

# A rate's SE from the mean of the two estimates, kept off 0 and 1 (two rates of exactly 0 or 1 stay equal).
rate_se <- function(column, found, reference) {
  rate <- pmin(pmax((found + reference) / 2, 0.5 / iterations), 1 - 0.5 / iterations)
  sqrt(rate * (1 - rate) * 2 / iterations)
}

replicate_kurtosis <- function(draws) {
  apply(draws, 2L, \(v) {
    v <- v[!is.na(v)]
    centred <- v - mean(v)
    mean(centred^4) / mean(centred^2)^2
  })
}

test_that("the bootstrap means and rates of every pathway equal the definitions on boot::boot() resamples (seeded Monte Carlo)", {
  skip_without_reference("boot")
  invisible(lapply(names(cases), \(name) {
    case <- cases[[name]]
    reference <- case$boot
    means <- standardised_pair(case$found, reference, paste0("mean_", statistics),
                               \(column, found, theirs) reference[[sub("^mean_", "sd_", column)]] * sqrt(2 / iterations))
    expect_reference(means$found, means$reference,
                     reference = boot_reference, covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: bootstrap means in standard errors of the difference", name), tolerance = 4,
                     key = "pathway",
                     note = "both sides divided by sd * sqrt(1/2000 + 1/2000), sd from the reference's replicates; 4 SE fixed before the run")
    rates <- standardised_pair(case$found, reference, rate_columns, rate_se)
    expect_reference(rates$found, rates$reference,
                     reference = boot_reference, covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: stability, informative, flip-consistency rates and p_stability in standard errors", name),
                     tolerance = 4, key = "pathway",
                     note = "both sides divided by sqrt(r (1 - r) (1/2000 + 1/2000)), r the mean of the two; 4 SE fixed before the run")
  }))
})

test_that("the bootstrap SDs and percentile bounds of every pathway agree with the definitions on boot::boot() resamples (seeded Monte Carlo)", {
  skip_without_reference("boot")
  invisible(lapply(names(cases), \(name) {
    case <- cases[[name]]
    reference <- case$boot
    draws <- attr(reference, "replicates")
    kurtosis <- lapply(draws, replicate_kurtosis)
    sds <- standardised_pair(case$found, reference, paste0("sd_", statistics), \(column, found, theirs) {
      statistic <- sub("^sd_", "", column)
      theirs * sqrt((kurtosis[[statistic]] - 1) / 4 * 2 / iterations)
    })
    expect_reference(sds$found, sds$reference, reference = boot_reference,
                     covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: bootstrap SDs in standard errors of the difference", name), tolerance = 4,
                     key = "pathway",
                     note = "both sides divided by sd * sqrt((kurtosis - 1) / 4 * (1/2000 + 1/2000)), from the reference's replicates; 4 SE fixed before the run")
    rows <- match(reference$pathway, case$found$pathway)
    # The reference's CDF at the package's bound: the point of [F(x-), F(x)] nearest the nominal level.
    level_at <- function(bounds, statistic, level) {
      vapply(seq_along(bounds), \(j) {
        values <- draws[[statistic]][, j]
        values <- values[!is.na(values)]
        if (is.na(bounds[[j]]) || !length(values)) return(NA_real_)
        min(max(level, mean(values < bounds[[j]])), mean(values <= bounds[[j]]))
      }, numeric(1))
    }
    bounds <- expand.grid(statistic = statistics, side = c("lo", "hi"), stringsAsFactors = FALSE)
    levels <- do.call(rbind, lapply(seq_len(nrow(bounds)), \(i) {
      statistic <- bounds$statistic[[i]]
      level <- if (bounds$side[[i]] == "lo") 0.025 else 0.975
      se <- sqrt(level * (1 - level) * 2 / iterations)
      found <- case$found[[sprintf("ci_%s_%s", statistic, bounds$side[[i]])]][rows]
      # Divergence and G2 are not defined at the root, which has no parent.
      defined <- !(statistic %in% c("divergence", "G2") & reference$pathway == "(start)")
      data.frame(pathway = reference$pathway, bound = sprintf("ci_%s_%s", statistic, bounds$side[[i]]),
                 found = level_at(found, statistic, level) / se, expected = ifelse(defined, level / se, NA_real_),
                 stringsAsFactors = FALSE)
    }))
    expect_reference(data.frame(pathway = levels$pathway, bound = levels$bound, level = levels$found),
                     data.frame(pathway = levels$pathway, bound = levels$bound, level = levels$expected),
                     reference = boot_reference, covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: reference CDF at the package's 95%% percentile bounds, in standard errors", name),
                     tolerance = 4, key = c("pathway", "bound"),
                     note = "level of the bound in the reference's 2000 replicates against 0.025 / 0.975, both divided by sqrt(a (1 - a) (1/2000 + 1/2000)); 4 SE fixed before the run")
  }))
})

test_that("as.data.frame() of the bootstrap is its summary, the same rates as the definitions on boot::boot() resamples", {
  skip_without_reference("boot")
  case <- cases[[1L]]
  found <- as.data.frame(bootstrap_pathways(context_tree(case$sequences, max_depth = 3L, min_count = 5L,
                                                         alphabet = alphabet), iter = iterations, seed = 1L))
  rates <- standardised_pair(found, case$boot, rate_columns, rate_se)
  expect_reference(rates$found, rates$reference,
                   reference = boot_reference, covers = c("bootstrap_pathways()", "as.data.frame.transitiontrees_bootstrap()"),
                   quantity = "seed 1: rates from as.data.frame() in standard errors", tolerance = 4, key = "pathway",
                   note = "both sides divided by sqrt(r (1 - r) (1/2000 + 1/2000)), r the mean of the two; 4 SE fixed before the run")
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
                 G2 = if (pathway == "(start)") NA_real_ else {
                   counts <- as.numeric(all_counts[pathway, ])
                   q <- distribution(ref_parent(pathway))
                   seen <- counts > 0
                   2 * sum(counts[seen] * log(counts[seen] / (sum(counts) * q[seen])))
                 },
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
    expect_reference(case$found, subset(expected, select = c(pathway, G2)),
                     reference = "independent: G2 of PST::cprob() counts against the parent's distribution",
                     covers = c("bootstrap_pathways()", "summary.transitiontrees_bootstrap()"),
                     quantity = sprintf("%s: observed G2 against the parent", name), tolerance = 1e-10,
                     measure = "relative", key = "pathway",
                     note = "2 sum n_s log(n_s / (n q_s)) over the states seen, PST::cprob() counts of the pathway and its parent")
  }))
})
