# Bounded generation supplements deterministic regressions. Finite pair sets
# determine a paired mean/SE/t interval independently of row and ID labels.
# Three complete profiles guarantee estimability; chronology and detection are
# fixed. BH is calculated from sorted raw p-values, without production helpers.
test_that("bounded finite-pair cases match arithmetic and preserve subject pairing", {
  withr::local_preserve_seed()
  withr::local_seed(20261012)
  for (case in seq_len(12L)) {
    n <- sample(6:12, 1L)
    k <- sample(2:5, 1L)
    values <- matrix(round(stats::rnorm(n * k, 20, 4), 3), n, k)
    unavailable <- sample(seq.int(4L, n * k), max(1L, floor(n * k / 6)))
    # Only incomplete subjects are perturbed, preserving three whole profiles.
    unavailable <- unavailable[((unavailable - 1L) %% n) >= 3L]
    values[unavailable] <- rep(c(NA_real_, NaN, Inf, -Inf), length.out = length(unavailable))
    vars <- paste0("score_t", seq_len(k) - 1L)
    data <- data.frame(subject_id = paste0("S", seq_len(n)), values,
                       check.names = FALSE)
    names(data)[-1L] <- vars
    args <- list(id = "subject_id", outcomes = "score", time_vars = vars,
                 covariates = character(), analyses = "correlations",
                 alpha = 0.1, p_adjust_method = "BH", verbose = FALSE)
    label <- sprintf("seed=20261012 case=%d n=%d k=%d", case, n, k)
    non_finite <- sum(!is.na(values) & !is.finite(values))
    warning_message <- paste0(non_finite,
      " Inf/-Inf values in longitudinal variables will be treated as NA in analyses.")
    if (non_finite > 0L) {
      expect_warning(result <- do.call(mira_test_info_base, c(list(data = data), args)),
                     warning_message, fixed = TRUE)
    } else {
      expect_no_warning(result <- do.call(mira_test_info_base, c(list(data = data), args)))
    }
    finite <- is.finite(values)
    expect_equal(result$descriptives$n, colSums(finite), info = label)
    expect_equal(unname(result$correlations$pairwise_n), unname(crossprod(finite)), info = label)
    pairs <- utils::combn(seq_len(k), 2L)
    expected_p <- numeric(ncol(pairs))
    for (i in seq_len(ncol(pairs))) {
      a <- pairs[1L, i]
      b <- pairs[2L, i]
      keep <- finite[, a] & finite[, b]
      delta <- values[keep, b] - values[keep, a]
      count <- length(delta)
      center <- sum(delta) / count
      sd_delta <- sqrt(sum((delta - center)^2) / (count - 1L))
      se_delta <- sd_delta / sqrt(count)
      half_width <- stats::qt(0.95, count - 1L) * se_delta
      expected_p[i] <- 2 * stats::pt(-abs(center / se_delta), count - 1L)
      observed <- result$change[i, ]
      expect_identical(c(observed$from, observed$to), vars[c(a, b)], info = label)
      expect_equal(observed$n, count, info = label)
      expect_equal(c(observed$mean_change, observed$sd_change, observed$se_change),
                   c(center, sd_delta, se_delta), tolerance = 1e-12, info = label)
      expect_equal(c(observed$ci_lower, observed$ci_upper),
                   center + c(-half_width, half_width), tolerance = 1e-12, info = label)
      expect_equal(observed$paired_t_p, expected_p[i], tolerance = 1e-12, info = label)
    }
    order_p <- order(expected_p)
    adjusted <- numeric(length(expected_p))
    adjusted[order_p] <- pmin(1, rev(cummin(rev(
      expected_p[order_p] * length(expected_p) / seq_along(expected_p)))))
    expect_equal(result$change$paired_t_p_adj, adjusted, tolerance = 1e-12, info = label)
    # ID is a key, not a numeric regressor. Renaming it bijectively and
    # permuting whole subject rows leaves every paired contrast unchanged.
    permutation <- sample(seq_len(n))
    renamed <- data[permutation, c("subject_id", rev(vars)), drop = FALSE]
    renamed$subject_id <- paste0("renamed_", permutation * 17L)
    if (non_finite > 0L) {
      expect_warning(changed <- do.call(mira_test_info_base, c(list(data = renamed), args)),
                     warning_message, fixed = TRUE)
    } else {
      expect_no_warning(changed <- do.call(mira_test_info_base, c(list(data = renamed), args)))
    }
    expect_equal(changed$change, result$change, tolerance = 1e-12, info = label)
    expect_equal(changed$descriptives, result$descriptives, tolerance = 1e-12, info = label)
  }
})
