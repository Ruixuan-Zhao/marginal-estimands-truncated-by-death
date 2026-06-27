# ===============================================================
# Simulation wrapper for while guaranteed-survival and while 
# extended-survival estimators
# ===============================================================

#' Compute truth values of \eqn{\mu(1) - \mu(0)} and \eqn{\mu^{\text{ext}}}
#'
#' @keywords internal
#' @noRd
make_truth_vectors <- function(tau_time = c(0, 1/4, 1/2, 1),
                               n.truth = 1000000,
                               seed = 12345) {
  truth_dat <- simulate_SubA(
    n.sample = n.truth,
    seed = seed,
    keep_potential = TRUE
  )
  
  truth_obj <- compute_truth_values(
    potential_data = truth_dat$potential_data,
    tau_time = tau_time
  )
  
  true_diffmu <- setNames(
    as.numeric(truth_obj$mu_diff$mu_diff),
    truth_obj$mu_diff$weight_type
  )
  
  true_mu.ext <- setNames(
    as.numeric(truth_obj$mu_ext$mu_ext),
    truth_obj$mu_ext$phi_type
  )
  
  list(
    diffmu = true_diffmu,
    mu.ext = true_mu.ext
  )
}


#' Summarize bias and Monte Carlo standard error
#'
#' @keywords internal
#' @noRd
summarize_estimates <- function(est_mat,
                                truth_vec,
                                n.sample) {
  rows <- lapply(colnames(est_mat), function(type) {
    est <- est_mat[, type]
    est <- est[is.finite(est)]
    truth <- as.numeric(truth_vec[type])
    
    if (length(est) == 0 || !is.finite(truth)) {
      bias <- NA_real_
      bias_se <- NA_real_
    } else {
      bias_vec <- est - truth
      bias <- mean(bias_vec)
      bias_se <- sd(bias_vec) / sqrt(length(bias_vec))
    }
    
    data.frame(
      type = type,
      n.sample = n.sample,
      true_value = truth,
      bias = bias,
      bias_se = bias_se,
      row.names = NULL
    )
  })
  
  do.call(rbind, rows)
}


#' Run the simulation for one sample size
#'
#' @keywords internal
#' @noRd
run_sim_SubA_one_setting <- function(n.sample,
                                     MC,
                                     truth,
                                     tau_time = c(0, 1/4, 1/2, 1),
                                     seed = 12345,
                                     verbose = TRUE) {
  diff_names <- names(truth$diffmu)
  mu.ext_names <- names(truth$mu.ext)
  
  est_diffmu_mat <- matrix(
    NA_real_,
    nrow = MC,
    ncol = length(diff_names),
    dimnames = list(NULL, diff_names)
  )
  
  est_mu.ext_mat <- matrix(
    NA_real_,
    nrow = MC,
    ncol = length(mu.ext_names),
    dimnames = list(NULL, mu.ext_names)
  )
  
  set.seed(seed)
  sim_seeds <- sample.int(1e8, MC)
  
  for (mc in seq_len(MC)) {
    if (verbose) {
      print_now <- (mc == 1) || (mc %% max(1, floor(MC / 10)) == 0)
      if (print_now) {
        cat("n =", n.sample, ", MC =", mc, "/", MC, "\n")
      }
    }
    
    est <- tryCatch(
      {
        dat <- simulate_SubA(
          n.sample = n.sample,
          seed = sim_seeds[mc]
        )
        
        get_est_SubA_allWeight(
          data = dat$data,
          tau_time = tau_time
        )
      },
      error = function(e) NULL
    )
    
    if (is.null(est)) {
      next
    }
    
    est_diffmu_mat[mc, diff_names] <- est$Diff_mu[diff_names]
    est_mu.ext_mat[mc, mu.ext_names] <- est$mu.ext[mu.ext_names]
  }
  
  list(
    summary_diffmu = summarize_estimates(
      est_mat = est_diffmu_mat,
      truth_vec = truth$diffmu,
      n.sample = n.sample
    ),
    summary_mu.ext = summarize_estimates(
      est_mat = est_mu.ext_mat,
      truth_vec = truth$mu.ext,
      n.sample = n.sample
    )
  )
}


#' Run the simulation over a grid of sample sizes
#'
#' Computes simulation summaries for the while guaranteed-survival contrast
#' \eqn{\mu(1) - \mu(0)} and the while extended-survival estimand
#' \eqn{\mu^{\text{ext}}}. For each sample size, the function runs `MC`
#' simulated data sets, compares estimates with truth values, and 
#' returns bias and Monte Carlo standard error summaries.
#'
#' @param n.sample.vec Numeric vector of sample sizes.
#' @param MC Number of Monte Carlo repetitions for each sample size.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#' @param seed Integer random seed used to generate truth and simulation seeds.
#' @param n.truth Sample size used to approximate the truth values.
#' @param verbose Logical. If `TRUE`, print simulation progress.
#'
#' @return A list with:
#'   \describe{
#'     \item{`summary_diffmu`}{A data frame with simulation summaries for
#'       \eqn{\mu(1) - \mu(0)}.}
#'     \item{`summary_mu.ext`}{A data frame with simulation summaries for
#'       \eqn{\mu^{\text{ext}}}.}
#'   }
#'   Each summary data frame contains `type`, `n.sample`, `true_value`, `bias`,
#'   and `bias_se`.
#'
#' @export
run_sim_SubA_grid <- function(n.sample.vec,
                              MC,
                              tau_time = c(0, 1/4, 1/2, 1),
                              seed = 12345,
                              n.truth = 1000000,
                              verbose = TRUE) {
  if (verbose) {
    cat("Computing truth with n.truth =", n.truth, "...\n")
  }
  
  truth <- make_truth_vectors(
    tau_time = tau_time,
    n.truth = n.truth,
    seed = seed + 999
  )
  
  summary_diffmu_all <- data.frame()
  summary_mu.ext_all <- data.frame()
  
  for (i in seq_along(n.sample.vec)) {
    n_i <- n.sample.vec[i]
    MC_i <- MC
    
    if (verbose) {
      cat("\n============================================================\n")
      cat("Running simulation: n =", n_i, ", MC =", MC_i, "\n")
      cat("============================================================\n")
    }
    
    res_i <- run_sim_SubA_one_setting(
      n.sample = n_i,
      MC = MC_i,
      truth = truth,
      tau_time = tau_time,
      seed = seed + 1000 * i,
      verbose = verbose
    )
    
    summary_diffmu_all <- rbind(summary_diffmu_all, res_i$summary_diffmu)
    summary_mu.ext_all <- rbind(summary_mu.ext_all, res_i$summary_mu.ext)
  }
  
  rownames(summary_diffmu_all) <- NULL
  rownames(summary_mu.ext_all) <- NULL
  
  list(
    summary_diffmu = summary_diffmu_all,
    summary_mu.ext = summary_mu.ext_all
  )
}
