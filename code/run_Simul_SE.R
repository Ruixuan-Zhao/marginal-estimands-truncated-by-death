#' Run one simulation setting for marginal separable effect estimands
#'
#' @keywords internal
#' @noRd
run_sim_SE_one_setting <- function(n.sample,
                                   MC,
                                   truth,
                                   tau_time = c(0, 1 / 4, 1 / 2, 1),
                                   seed = 12345,
                                   verbose = TRUE) {
  set.seed(seed)

  weight_names <- rownames(truth)
  contrast_names <- colnames(truth)

  est_diff_gamma <- array(
    NA_real_,
    dim = c(MC, length(weight_names), length(contrast_names)),
    dimnames = list(NULL, weight_names, contrast_names)
  )

  sim_seeds <- sample.int(1e8, MC)

  for (mc in seq_len(MC)) {
    if (verbose && (mc %% max(1, floor(MC / 10)) == 0 || mc == 1)) {
      cat("n =", n.sample, ", MC =", mc, "/", MC, "\n")
    }

    fit_mc <- tryCatch({
      dat <- simulate_SE(
        n.sample = n.sample,
        seed = sim_seeds[mc]
      )

      get_est_Separable_Effect_allWeight(
        data = dat,
        tau_time = tau_time
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(fit_mc)) {
      next
    }

    est_diff_gamma[mc, , ] <-
      fit_mc$Diff_Gamma[weight_names, contrast_names, drop = FALSE]
  }

  summary_rows <- list()
  row_id <- 1

  for (weight_type in weight_names) {
    for (contrast in contrast_names) {
      est <- est_diff_gamma[, weight_type, contrast]
      est <- est[is.finite(est)]
      truth_value <- truth[weight_type, contrast]

      summary_rows[[row_id]] <- data.frame(
        type = weight_type,
        contrast = contrast,
        n.sample = n.sample,
        truth = as.numeric(truth_value),
        bias = mean(est - truth_value),
        bias_se = sd(est - truth_value) / sqrt(length(est)),
        row.names = NULL
      )
      row_id <- row_id + 1
    }
  }

  summary_diff_gamma <- do.call(rbind, summary_rows)

  list(summary_diff_gamma = summary_diff_gamma)
}


#' Run the marginal separable effect simulation over sample sizes
#'
#' @param n.sample.vec Numeric vector of sample sizes.
#' @param MC Number of Monte Carlo repetitions for each sample size.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#' @param seed Integer random seed used to generate truth and simulation seeds.
#' @param n.truth Sample size used to approximate the truth.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return A list with `summary_diff_gamma`, a data frame containing the
#'   columns `type`, `contrast`, `n.sample`, `truth`, `bias`, and `bias_se`.
#'
#' @export
run_sim_SE_grid <- function(n.sample.vec,
                            MC,
                            tau_time = c(0, 1 / 4, 1 / 2, 1),
                            seed = 12345,
                            n.truth = 1000000,
                            verbose = TRUE) {
  if (verbose) {
    cat("Computing SE truth with n.truth =", n.truth, "...\n")
  }

  truth <- compute_truth_values_SE(
    n.sample = n.truth,
    seed = seed + 999
  )$trueValue_SE

  summary_diff_gamma_all <- data.frame()

  for (i in seq_along(n.sample.vec)) {
    n_i <- n.sample.vec[i]

    if (verbose) {
      cat("\n============================================================\n")
      cat("Running SE simulation: n =", n_i, ", MC =", MC, "\n")
      cat("============================================================\n")
    }

    res_i <- run_sim_SE_one_setting(
      n.sample = n_i,
      MC = MC,
      truth = truth,
      tau_time = tau_time,
      seed = seed + 1000 * i,
      verbose = verbose
    )

    summary_diff_gamma_all <- rbind(summary_diff_gamma_all, res_i$summary_diff_gamma)
  }

  list(summary_diff_gamma = summary_diff_gamma_all)
}
