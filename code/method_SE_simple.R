################################################################################
###### Estimation method for marginal separable effect estimands and CSEs ######
################################################################################

#' Estimate marginal separable effect estimands
#'
#' Estimates marginal separable effect contrasts
#' \eqn{\Gamma(1, z_S) - \Gamma(0, z_S)}, the contrast
#' \eqn{\Gamma(1, 1) - \Gamma(1, 0)}, and the conditional separable effects
#' CSE(t), for \eqn{t = 1, \ldots, t_{\max}}, from a longitudinal data set.
#'
#' @param data A data frame with one row per subject. Required columns are
#'   `Z`, the treatment indicator; `SURVT`, the observed exit time; `Y0`, the
#'   baseline outcome; and time-indexed columns `S1`-`S[t_max]`,
#'   `L1`-`L[t_max]`, and `Y1`-`Y[t_max]`, where `t_max` is the number of
#'   follow-up times. All other columns are treated as baseline covariates.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#' @param use_L_history Logical. If `TRUE`, use progression history through
#'   time `r` in the outcome regression; if `FALSE`, use only progression at
#'   time `r`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`Diff_Gamma`}{Estimates of
#'       \eqn{\Gamma(1, z_S) - \Gamma(0, z_S)} under the exit-time, average,
#'       cumulative, and AUC-based weights.}
#'     \item{`Diff_Zs_given_Zy1`}{Estimates of \eqn{\Gamma(1, 1) - \Gamma(1, 0)}.}
#'     \item{`Diff_CSE`}{Estimates of CSE(t), for \eqn{t = 1, \ldots, t_{\max}}.}
#'   }
#'
#' @export
get_est_Separable_Effect_allWeight <- function(data,
                                                tau_time,
                                                use_L_history = TRUE) {
  t_max <- length(tau_time) - 1
  s_cols <- paste0("S", seq_len(t_max))
  y_cols <- paste0("Y", seq_len(t_max))
  l_cols <- paste0("L", seq_len(t_max))
  special_cols <- c("Z", "SURVT", "Y0", s_cols, l_cols, y_cols)
  baseline_covariate_cols <- setdiff(names(data), special_cols)

  Z <- data$Z
  L0 <- as.matrix(data[baseline_covariate_cols])
  Y0 <- data$Y0
  SURVT <- data$SURVT

  # Treatment model Pr(Z = z | L0).
  design_z <- cbind(Intercept = 1, L0)
  fit_z <- glm.fit(
    x = design_z,
    y = Z,
    family = binomial()
  )
  theta_z <- fit_z$coefficients
  p_Z1 <- plogis(drop(design_z %*% theta_z))
  p_Z0 <- 1 - p_Z1

  survival_prob <- compute_St_Zz(
    data = data,
    t_max = t_max,
    covariate_cols = baseline_covariate_cols
  )

  weight.sum <- Weights_general(tau_time)
  weights_time <- Weights_timeSpecified(
    survival_prob$Prob_St_Z0,
    survival_prob$Prob_St_Z1
  )
  weights_CSE_z0 <- weights_time$weights_timeSpecified_z0
  weights_CSE_z1 <- weights_time$weights_timeSpecified_z1

  Gamma_zy_zs <- matrix(0, 2, 2)
  colnames(Gamma_zy_zs) <- c("z_s=0", "z_s=1")
  rownames(Gamma_zy_zs) <- c("z_y=0", "z_y=1")

  Gamma_zy_zs.sum <- list(
    Gamma_Exit = Gamma_zy_zs,
    Gamma_Ave = Gamma_zy_zs,
    Gamma_Cum = Gamma_zy_zs,
    Gamma_AUC = Gamma_zy_zs
  )
  names(Gamma_zy_zs.sum) <- c("Exit time", "Average", "Cumulative", "AUC-based")

  Gamma_zy_zs.SACE <- setNames(
    replicate(t_max, Gamma_zy_zs, simplify = FALSE),
    paste("CSE at time", seq_len(t_max))
  )

  # mu_ab_T_b stores E[I(T(b) = t) Y^r(a, b)] for 0 <= r <= t <= t_max.
  mu_00_T_0 <- .make_components(t_max)
  mu_10_T_0 <- .make_components(t_max)
  mu_01_T_1 <- .make_components(t_max)
  mu_11_T_1 <- .make_components(t_max)

  for (t in 0:t_max) {
    mu_00_T_0[[t + 1]][1] <- mean((Z == 0) * (SURVT == t) * Y0 / p_Z0)
    mu_10_T_0[[t + 1]][1] <- mean((Z == 0) * (SURVT == t) * Y0 / p_Z0)
    mu_01_T_1[[t + 1]][1] <- mean((Z == 1) * (SURVT == t) * Y0 / p_Z1)
    mu_11_T_1[[t + 1]][1] <- mean((Z == 1) * (SURVT == t) * Y0 / p_Z1)
  }

  for (r in 1:t_max) {
    y <- data[[y_cols[r]]]
    S <- data[[s_cols[r]]]
    selected_l_cols <- if (use_L_history) l_cols[1:r] else l_cols[r]
    L <- as.matrix(data[, selected_l_cols, drop = FALSE])

    # Outcome model E(Y^r | S^r = 1, Z, L0, L).
    design_y <- cbind(
      Intercept = 1,
      L0,
      L,
      Z = Z
    )
    survivor_idx <- S == 1
    
    fit_y <- lm.fit(
      x = design_y[survivor_idx, , drop = FALSE],
      y = y[survivor_idx]
    )
    beta_y <- fit_y$coefficients
    beta_y[is.na(beta_y)] <- 0

    design_y_z0 <- cbind(
      Intercept = 1,
      L0,
      L,
      Z = 0
    )
    design_y_z1 <- cbind(
      Intercept = 1,
      L0,
      L,
      Z = 1
    )
    
    mu_y_z0 <- drop(design_y_z0 %*% beta_y)
    mu_y_z1 <- drop(design_y_z1 %*% beta_y)

    for (t in r:t_max) {
      mu_00_T_0[[t + 1]][r + 1] <- mean((Z == 0) * (SURVT == t) * mu_y_z0 / p_Z0)
      mu_10_T_0[[t + 1]][r + 1] <- mean((Z == 0) * (SURVT == t) * mu_y_z1 / p_Z0)
      mu_01_T_1[[t + 1]][r + 1] <- mean((Z == 1) * (SURVT == t) * mu_y_z0 / p_Z1)
      mu_11_T_1[[t + 1]][r + 1] <- mean((Z == 1) * (SURVT == t) * mu_y_z1 / p_Z1)
    }
  }

  for (m in 1:4) {
    for (t in 0:t_max) {
      Gamma_zy_zs.sum[[m]][1, 1] <- Gamma_zy_zs.sum[[m]][1, 1] +
        sum(weight.sum[[m]][[t + 1]] * mu_00_T_0[[t + 1]])
      Gamma_zy_zs.sum[[m]][1, 2] <- Gamma_zy_zs.sum[[m]][1, 2] +
        sum(weight.sum[[m]][[t + 1]] * mu_01_T_1[[t + 1]])
      Gamma_zy_zs.sum[[m]][2, 1] <- Gamma_zy_zs.sum[[m]][2, 1] +
        sum(weight.sum[[m]][[t + 1]] * mu_10_T_0[[t + 1]])
      Gamma_zy_zs.sum[[m]][2, 2] <- Gamma_zy_zs.sum[[m]][2, 2] +
        sum(weight.sum[[m]][[t + 1]] * mu_11_T_1[[t + 1]])
    }
  }

  for (h in 1:t_max) {
    for (t in 0:t_max) {
      Gamma_zy_zs.SACE[[h]][1, 1] <- Gamma_zy_zs.SACE[[h]][1, 1] +
        sum(weights_CSE_z0[[h]][[t + 1]] * mu_00_T_0[[t + 1]])
      Gamma_zy_zs.SACE[[h]][1, 2] <- Gamma_zy_zs.SACE[[h]][1, 2] +
        sum(weights_CSE_z1[[h]][[t + 1]] * mu_01_T_1[[t + 1]])
      Gamma_zy_zs.SACE[[h]][2, 1] <- Gamma_zy_zs.SACE[[h]][2, 1] +
        sum(weights_CSE_z0[[h]][[t + 1]] * mu_10_T_0[[t + 1]])
      Gamma_zy_zs.SACE[[h]][2, 2] <- Gamma_zy_zs.SACE[[h]][2, 2] +
        sum(weights_CSE_z1[[h]][[t + 1]] * mu_11_T_1[[t + 1]])
    }
  }

  Diff_Gamma <- matrix(0, 4, 2)
  colnames(Diff_Gamma) <- c(
    "Gamma(1,0)-Gamma(0,0)(z_s=0)",
    "Gamma(1,1)-Gamma(0,1)(z_s=1)"
  )
  rownames(Diff_Gamma) <- c("Exit time", "Average", "Cumulative", "AUC-based")

  for (m in 1:4) {
    Diff_Gamma[m, ] <- c(
      Gamma_zy_zs.sum[[m]][2, 1] - Gamma_zy_zs.sum[[m]][1, 1],
      Gamma_zy_zs.sum[[m]][2, 2] - Gamma_zy_zs.sum[[m]][1, 2]
    )
  }

  Diff_Zs_given_Zy1 <- c(
    "Cumulative" = Gamma_zy_zs.sum[["Cumulative"]][2, 2] -
      Gamma_zy_zs.sum[["Cumulative"]][2, 1],
    "AUC-based" = Gamma_zy_zs.sum[["AUC-based"]][2, 2] -
      Gamma_zy_zs.sum[["AUC-based"]][2, 1]
  )

  Diff_CSE <- matrix(0, t_max, 2)
  colnames(Diff_CSE) <- colnames(Diff_Gamma)
  rownames(Diff_CSE) <- paste("CSE at time", 1:t_max)

  for (m in 1:t_max) {
    Diff_CSE[m, ] <- c(
      Gamma_zy_zs.SACE[[m]][2, 1] - Gamma_zy_zs.SACE[[m]][1, 1],
      Gamma_zy_zs.SACE[[m]][2, 2] - Gamma_zy_zs.SACE[[m]][1, 2]
    )
  }

  list(
    Diff_Gamma = Diff_Gamma,
    Diff_Zs_given_Zy1 = Diff_Zs_given_Zy1,
    Diff_CSE = Diff_CSE
  )
}
