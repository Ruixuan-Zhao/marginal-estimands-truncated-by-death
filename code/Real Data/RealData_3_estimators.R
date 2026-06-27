################################################################################
######################## Real data analysis functions ##########################
################################################################################

correct_mu_ext <- function(original_mu.ext, data, tau_time) {
  t_max <- length(tau_time) - 1
  mu.ext_names <- names(original_mu.ext)
  if (is.null(mu.ext_names)) {
    mu.ext_names <- c("Cumulative", "AUC-based")[seq_along(original_mu.ext)]
  }

  non_covariate_cols <- c(
    "Z", "SURVT",
    paste0("Y", 0:t_max),
    paste0("S", seq_len(t_max)),
    paste0("L", seq_len(t_max))
  )
  covariate_cols <- setdiff(names(data), non_covariate_cols)

  Z <- data$Z
  A <- data$A
  SURVT <- data$SURVT

  tx <- .fit_treatment_prob(
    data = data,
    covariate_cols = covariate_cols
  )
  p_Z1 <- pmin(pmax(tx$p_Z1, 1e-8), 1 - 1e-8)
  p_Z0 <- pmin(pmax(tx$p_Z0, 1e-8), 1 - 1e-8)

  mu_A_T1 <- rep(0, t_max + 1)
  mu_A_T0 <- rep(0, t_max + 1)
  Prob_T_1 <- rep(0, t_max + 1)
  Prob_T_0 <- rep(0, t_max + 1)

  for (tt in 0:t_max) {
    mu_A_T1[tt + 1] <- mean((Z == 1) * (SURVT == tt) * A / p_Z1)
    mu_A_T0[tt + 1] <- mean((Z == 0) * (SURVT == tt) * A / p_Z0)
    Prob_T_1[tt + 1] <- mean((Z == 1) * (SURVT == tt) / p_Z1)
    Prob_T_0[tt + 1] <- mean((Z == 0) * (SURVT == tt) / p_Z0)
  }

  weight.sum <- Weights_general(tau_time)
  correction <- setNames(rep(0, length(mu.ext_names)), mu.ext_names)

  for (weight_name in mu.ext_names) {
    for (tt in 0:t_max) {
      mean_diff <- mu_A_T1[tt + 1] - mu_A_T0[tt + 1]
      correction[weight_name] <- correction[weight_name] +
        sum(weight.sum[[weight_name]][[tt + 1]]) * mean_diff
    }
  }

  corrected_mu.ext <- original_mu.ext + correction[mu.ext_names]
  names(corrected_mu.ext) <- mu.ext_names

  list(
    corrected_mu.ext = corrected_mu.ext,
    ExtendedTime = sum((Prob_T_1 - Prob_T_0) * tau_time)
  )
}


PointEst_all <- function(data, tau_time, use_L_history = FALSE) {
  t_max <- length(tau_time) - 1

  whilealive_data <- data[
    setdiff(names(data), c(paste0("S", seq_len(t_max)), paste0("L", seq_len(t_max))))
  ]
  re_WhileAlive <- get_est_WhileAlive_allWeight(
    data = whilealive_data,
    tau_time = tau_time
  )
  results.WhileAlive <- re_WhileAlive$Diff_lambda

  re_SubA <- get_est_SubA_allWeight(
    data = data,
    tau_time = tau_time,
    use_L_history = use_L_history
  )
  results.mu <- re_SubA$Diff_mu
  results.SACE <- estimate_SACE(data = data)
  results.mu.ext <- re_SubA$mu.ext
  corrected_mu.ext <- correct_mu_ext(
    original_mu.ext = results.mu.ext,
    data = data,
    tau_time = tau_time
  )
  results.corrected_mu.ext <- corrected_mu.ext$corrected_mu.ext
  results.ExtendedTime <- corrected_mu.ext$ExtendedTime

  re_SE <- get_est_Separable_Effect_allWeight(
    data = data,
    tau_time = tau_time,
    use_L_history = use_L_history
  )
  results.SE_zs0 <- re_SE$Diff_Gamma[, 1]
  results.CSE <- re_SE$Diff_CSE[, 1]
  results.SE_mu.ext <- re_SE$Diff_Zs_given_Zy1
  results.SE_corrected_mu.ext <- correct_mu_ext(
    original_mu.ext = results.SE_mu.ext,
    data = data,
    tau_time = tau_time
  )$corrected_mu.ext

  list(
    results.SACE = results.SACE,
    results.CSE = results.CSE,
    results.WhileAlive = results.WhileAlive,
    results.mu = results.mu,
    results.mu.ext = results.mu.ext,
    results.corrected_mu.ext = results.corrected_mu.ext,
    results.ExtendedTime = results.ExtendedTime,
    results.SE_zs0 = results.SE_zs0,
    results.SE_mu.ext = results.SE_mu.ext,
    results.SE_corrected_mu.ext = results.SE_corrected_mu.ext
  )
}


Bootstrapping <- function(data, tau_time, rep.each = 500, use_L_history = FALSE) {
  t_max <- length(tau_time) - 1

  results.SACE <- matrix(NA, rep.each, t_max)
  colnames(results.SACE) <- paste("SACE at time", 1:t_max)
  results.CSE <- matrix(NA, rep.each, t_max)
  colnames(results.CSE) <- paste("CSE at time", 1:t_max)

  results.WhileAlive <- matrix(NA, rep.each, 4)
  colnames(results.WhileAlive) <- c("Exit time", "Average", "Cumulative", "AUC-based")
  results.mu <- matrix(NA, rep.each, 4)
  colnames(results.mu) <- c("Exit time", "Average", "Cumulative", "AUC-based")
  results.mu.ext <- matrix(NA, rep.each, 2)
  colnames(results.mu.ext) <- c("Cumulative", "AUC-based")
  results.corrected_mu.ext <- matrix(NA, rep.each, 2)
  colnames(results.corrected_mu.ext) <- c("Cumulative", "AUC-based")
  results.ExtendedTime <- matrix(NA, rep.each, 1)
  colnames(results.ExtendedTime) <- "E[sum tau_t I(t=T(1))]-E[sum tau_t I(t=T(0))]"
  results.SE_zs0 <- matrix(NA, rep.each, 4)
  colnames(results.SE_zs0) <- c("Exit time", "Average", "Cumulative", "AUC-based")
  results.SE_mu.ext <- matrix(NA, rep.each, 2)
  colnames(results.SE_mu.ext) <- c("Cumulative", "AUC-based")
  results.SE_corrected_mu.ext <- matrix(NA, rep.each, 2)
  colnames(results.SE_corrected_mu.ext) <- c("Cumulative", "AUC-based")

  for (rep in seq_len(rep.each)) {
    message("Bootstrap replication ", rep)
    set.seed(rep)

    n <- nrow(data)
    index <- sample(seq_len(n), n, replace = TRUE)
    data.B <- data[index, , drop = FALSE]
    res_points <- PointEst_all(
      data = data.B,
      tau_time = tau_time,
      use_L_history = use_L_history
    )

    results.SACE[rep, ] <- res_points$results.SACE
    results.CSE[rep, ] <- res_points$results.CSE
    results.WhileAlive[rep, ] <- res_points$results.WhileAlive
    results.mu[rep, ] <- res_points$results.mu
    results.mu.ext[rep, ] <- res_points$results.mu.ext
    results.corrected_mu.ext[rep, ] <- res_points$results.corrected_mu.ext
    results.ExtendedTime[rep, ] <- res_points$results.ExtendedTime
    results.SE_zs0[rep, ] <- res_points$results.SE_zs0
    results.SE_mu.ext[rep, ] <- res_points$results.SE_mu.ext
    results.SE_corrected_mu.ext[rep, ] <- res_points$results.SE_corrected_mu.ext
  }

  list(
    results.SACE = results.SACE,
    results.CSE = results.CSE,
    results.WhileAlive = results.WhileAlive,
    results.mu = results.mu,
    results.mu.ext = results.mu.ext,
    results.corrected_mu.ext = results.corrected_mu.ext,
    results.ExtendedTime = results.ExtendedTime,
    results.SE_zs0 = results.SE_zs0,
    results.SE_mu.ext = results.SE_mu.ext,
    results.SE_corrected_mu.ext = results.SE_corrected_mu.ext
  )
}


get_boot_summary <- function(mat, point_est) {
  data.frame(
    Estimator = colnames(mat),
    Point_Estimate = point_est,
    Bootstrap_SE = apply(mat, 2, sd, na.rm = TRUE),
    CI_Lower = apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_Upper = apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    row.names = NULL
  )
}


run_realdata_analysis <- function(data,
                                  tau_time,
                                  boot_reps = 500,
                                  use_L_history = FALSE) {
  est_points <- PointEst_all(
    data = data,
    tau_time = tau_time,
    use_L_history = use_L_history
  )
  re_Boot <- Bootstrapping(
    data = data,
    tau_time = tau_time,
    rep.each = boot_reps,
    use_L_history = use_L_history
  )

  result_list <- lapply(names(re_Boot), function(nm) {
    get_boot_summary(re_Boot[[nm]], est_points[[nm]])
  })
  names(result_list) <- names(re_Boot)

  list(
    est_points = est_points,
    re_Boot = re_Boot,
    result_list = result_list
  )
}
