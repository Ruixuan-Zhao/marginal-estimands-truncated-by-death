# ============================================================
# Estimation methods for while guaranteed-survival and 
# while extended-survival estimands
# ============================================================

get_risk <- function(data, rr) {
  if (rr == 1) {
    rep(TRUE, nrow(data))
  } else {
    data[[paste0("S", rr - 1)]] == 1
  }
}


#' Negative log-likelihood for the survival model
#'
#' @keywords internal
#' @noRd
nLL_beta_gamma_homBG <- function(par, data, W, t_max, eps = 1e-8) {
  W_slope <- W[, -1, drop = FALSE]
  q <- ncol(W_slope)
  
  beta_int <- par[seq_len(t_max)]
  beta_slope <- par[t_max + seq_len(q)]
  
  offset <- t_max + q
  gamma_int <- par[offset + seq_len(t_max)]
  gamma_slope <- par[offset + t_max + seq_len(q)]
  
  loglike <- 0
  
  for (rr in seq_len(t_max)) {
    S <- as.numeric(data[[paste0("S", rr)]])
    Z <- as.numeric(data$Z)
    idx <- get_risk(data, rr)
    
    eta_beta <- beta_int[rr] + drop(W_slope %*% beta_slope)
    eta_gamma <- gamma_int[rr] + drop(W_slope %*% gamma_slope)
    
    p_z1 <- .clamp_prob(plogis(eta_beta), eps)
    p_z0 <- .clamp_prob(p_z1 * plogis(eta_gamma), eps)
    
    S_use <- S[idx]
    Z_use <- Z[idx]
    
    loglike <- loglike + sum(
      Z_use * (
        S_use * log(p_z1[idx]) +
          (1 - S_use) * log(1 - p_z1[idx])
      ) +
        (1 - Z_use) * (
          S_use * log(p_z0[idx]) +
            (1 - S_use) * log(1 - p_z0[idx])
        )
    )
  }
  
  nll <- -loglike
  if (!is.finite(nll)) nll <- 1e100
  
  nll
}


#' Construct starting values for the survival model
#'
#' @keywords internal
#' @noRd
make_start_beta_gamma <- function(data, W, t_max) {
  q <- ncol(W) - 1
  beta_int <- rep(0, t_max)
  gamma_int <- rep(0, t_max)
  
  for (rr in seq_len(t_max)) {
    S <- as.numeric(data[[paste0("S", rr)]])
    Z <- as.numeric(data$Z)
    risk <- get_risk(data, rr)
    
    p1_hat <- mean(S[risk & Z == 1])
    p0_hat <- mean(S[risk & Z == 0])
    
    p1_hat <- .clamp_prob(p1_hat, eps = 1e-4)
    p0_hat <- .clamp_prob(p0_hat, eps = 1e-4)
    
    beta_int[rr] <- qlogis(p1_hat)
    gamma_int[rr] <- qlogis(.clamp_prob(p0_hat / p1_hat, eps = 1e-4))
  }
  
  c(beta_int, rep(0, q), gamma_int, rep(0, q))
}


#' Estimate survival model parameters
#'
#' @keywords internal
#' @noRd
estimate_beta_gamma <- function(data, W, t_max) {
  q <- ncol(W) - 1
  
  opt <- optim(
    par = make_start_beta_gamma(data, W, t_max),
    fn = nLL_beta_gamma_homBG,
    data = data,
    W = W,
    t_max = t_max,
    method = "BFGS",
    control = list(maxit = 5000)
  )
  
  par_hat <- opt$par
  
  beta_int <- par_hat[seq_len(t_max)]
  beta_slope <- par_hat[t_max + seq_len(q)]
  
  offset <- t_max + q
  gamma_int <- par_hat[offset + seq_len(t_max)]
  gamma_slope <- par_hat[offset + t_max + seq_len(q)]
  
  beta <- vector("list", t_max)
  gamma <- vector("list", t_max)
  
  for (rr in seq_len(t_max)) {
    beta[[rr]] <- c(beta_int[rr], beta_slope)
    gamma[[rr]] <- c(gamma_int[rr], gamma_slope)
    
    names(beta[[rr]]) <- colnames(W)
    names(gamma[[rr]]) <- colnames(W)
  }
  
  list(beta = beta, gamma = gamma)
}


#' Fit the outcome regression at time r
#'
#' @keywords internal
#' @noRd
fit_alpha_r <- function(data,
                        rr,
                        gamma,
                        W,
                        baseline_covariate_cols,
                        use_L_history = TRUE) {
  y <- data[[paste0("Y", rr)]]
  Z <- data$Z
  l_cols <- if (use_L_history) {
    paste0("L", seq_len(rr))
  } else {
    paste0("L", rr)
  }
  Lbar <- as.matrix(data[l_cols])
  X <- as.matrix(data[baseline_covariate_cols])
  
  prod_gamma <- rep(1, nrow(data))
  
  for (k in seq_len(rr)) {
    prod_gamma <- prod_gamma * plogis(drop(W %*% gamma[[k]]))
  }
  
  X_alpha <- cbind(
    Intercept = 1,
    Lbar,
    X,
    A = data$A,
    Wexpit = (1 - prod_gamma) * Z,
    Z = Z
  )
  survivor_idx <- data[[paste0("S", rr)]] == 1
  
  fit <- lm.fit(
    x = X_alpha[survivor_idx, , drop = FALSE],
    y = y[survivor_idx]
  )
  
  alpha <- fit$coefficients
  alpha[is.na(alpha)] <- 0
  
  if (length(alpha) < ncol(X_alpha)) {
    alpha_full <- rep(0, ncol(X_alpha))
    alpha_full[seq_along(alpha)] <- alpha
    alpha <- alpha_full
  }
  
  X_alpha_LL1 <- cbind(
    Intercept = 1,
    Lbar,
    X,
    A = data$A,
    Wexpit = 0,
    Z = 1
  )
  
  mu_y_z1 <- rep(0, nrow(data))
  mu_y_z1[survivor_idx] <- drop(
    X_alpha_LL1[survivor_idx, , drop = FALSE] %*% alpha
  )
  
  mu_y_z1
}


#' Estimate the while guaranteed-survival and while extended-survival
#' estimands
#'
#' Estimates \eqn{\mu(1) - \mu(0)} and \eqn{\mu^{\text{ext}}} from a 
#' longitudinal data set. The function returns estimates for different weights.
#'
#' @param data A data frame with one row per subject. Required columns are
#'   `Z`, the treatment indicator; `A`, the substitution variable; `SURVT`,
#'   the observed exit time; `Y0`, the baseline outcome; baseline covariates
#'   named `X1`, `X2`, ..., `Xp`; and time-indexed columns `S1`-`S[t_max]`,
#'   `L1`-`L[t_max]`, and `Y1`-`Y[t_max]`, where `t_max` is the number of
#'   follow-up times.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#' @param use_L_history Logical. If `TRUE`, use progression history through
#'   time `r` in the outcome regression; if `FALSE`, use only progression at
#'   time `r`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`Diff_mu`}{Estimates of \eqn{\mu(1) - \mu(0)} under the
#'       exit-time, average, cumulative, and AUC-based weights.}
#'     \item{`mu.SACE`}{Estimates of SACE(t) for \eqn{t = 1, \ldots, t_{\max}}.}
#'     \item{`mu.ext`}{Estimates of \eqn{\mu^{\text{ext}}} under the
#'       cumulative and AUC-based weights.}
#'   }
#'
#' @export
get_est_SubA_allWeight <- function(data,
                                   tau_time = c(0, 1/4, 1/2, 1),
                                   use_L_history = TRUE) {
  t_max <- length(tau_time) - 1
  baseline_covariate_cols <- grep("^X[0-9]+$", names(data), value = TRUE)
  baseline_covariate_cols <- baseline_covariate_cols[
    order(as.integer(sub("^X", "", baseline_covariate_cols)))
  ]
  
  Z <- data$Z
  Y0 <- data$Y0
  SURVT <- data$SURVT
  W <- cbind(
    Intercept = 1,
    as.matrix(data[baseline_covariate_cols]),
    A = data$A
  )
  weight.sum <- Weights_general(tau_time)
  
  tx <- .fit_treatment_prob(
    data = data,
    covariate_cols = c(baseline_covariate_cols, "A")
  )
  
  p_Z1 <- .clamp_prob(tx$p_Z1)
  p_Z0 <- .clamp_prob(tx$p_Z0)
  
  survival_prob <- compute_St_Zz(
    data = data,
    t_max = t_max,
    covariate_cols = c("A", baseline_covariate_cols)
  )
  weights_naive <- Weights_timeSpecified(
    survival_prob$Prob_St_Z0,
    survival_prob$Prob_St_Z1
  )
  
  mu_0_T_0 <- vector("list", t_max + 1)
  mu_1_T_0 <- vector("list", t_max + 1)
  mu_ext_first <- vector("list", t_max + 1)
  
  for (tt in 0:t_max) {
    mu_0_T_0[[tt + 1]] <- rep(0, tt + 1)
    mu_1_T_0[[tt + 1]] <- rep(0, tt + 1)
    mu_ext_first[[tt + 1]] <- rep(0, tt + 1)
    
    mu_0_T_0[[tt + 1]][1] <- mean((Z == 0) * (SURVT == tt) * Y0 / p_Z0)
    mu_ext_first[[tt + 1]][1] <- mean((Z == 1) * (SURVT == tt) * Y0 / p_Z1)
  }
  
  beta_gamma <- estimate_beta_gamma(data, W, t_max)
  beta <- beta_gamma$beta
  gamma <- beta_gamma$gamma
  
  for (rr in seq_len(t_max)) {
    y_r <- data[[paste0("Y", rr)]]
    mu_y_z1 <- fit_alpha_r(
      data = data,
      rr = rr,
      gamma = gamma,
      W = W,
      baseline_covariate_cols = baseline_covariate_cols,
      use_L_history = use_L_history
    )
    
    for (tt in rr:t_max) {
      prod_gamma_t <- rep(1, nrow(data))
      
      for (k in seq_len(tt)) {
        prod_gamma_t <- prod_gamma_t * plogis(drop(W %*% gamma[[k]]))
      }
      
      if (tt == t_max) {
        pi_t <- prod_gamma_t
      } else {
        ebeta_next <- plogis(drop(W %*% beta[[tt + 1]]))
        egamma_next <- plogis(drop(W %*% gamma[[tt + 1]]))
        pi_t <- prod_gamma_t * (1 - ebeta_next * egamma_next)
      }
      
      obs_mu0 <- ifelse(
        is.na(y_r),
        0,
        (Z == 0) * (SURVT == tt) * y_r / p_Z0
      )
      
      mu_0_T_0[[tt + 1]][rr + 1] <- mean(obs_mu0)
      
      mu_1_T_0[[tt + 1]][rr + 1] <- mean(
        (Z == 1) * data[[paste0("S", tt)]] * pi_t * mu_y_z1 / p_Z1
      )
      
      obs_ext_first <- ifelse(
        is.na(y_r),
        0,
        (Z == 1) * (SURVT == tt) * y_r / p_Z1
      )
      
      mu_ext_first[[tt + 1]][rr + 1] <- mean(obs_ext_first)
    }
  }
  
  mu_01 <- lapply(weight.sum, function(x) c("mu(0)" = 0, "mu(1)" = 0))
  
  for (m in seq_along(weight.sum)) {
    for (tt in 0:t_max) {
      mu_01[[m]]["mu(0)"] <- mu_01[[m]]["mu(0)"] +
        sum(weight.sum[[m]][[tt + 1]] * mu_0_T_0[[tt + 1]])
      
      mu_01[[m]]["mu(1)"] <- mu_01[[m]]["mu(1)"] +
        sum(weight.sum[[m]][[tt + 1]] * mu_1_T_0[[tt + 1]])
    }
  }
  
  Diff_mu <- sapply(mu_01, function(x) x["mu(1)"] - x["mu(0)"])
  names(Diff_mu) <- names(weight.sum)
  
  mu.SACE <- rep(NA_real_, t_max)
  names(mu.SACE) <- paste("SACE at time", seq_len(t_max))
  
  for (rr in seq_len(t_max)) {
    weights_z0 <- weights_naive$weights_timeSpecified_z0[[rr]]
    mu0 <- 0
    mu1 <- 0
    
    for (tt in 0:t_max) {
      mu0 <- mu0 + sum(weights_z0[[tt + 1]] * mu_0_T_0[[tt + 1]])
      mu1 <- mu1 + sum(weights_z0[[tt + 1]] * mu_1_T_0[[tt + 1]])
    }
    
    mu.SACE[rr] <- mu1 - mu0
  }
  
  mu.ext <- c("Cumulative" = 0, "AUC-based" = 0)
  
  for (type in names(mu.ext)) {
    phi <- weight.sum[[type]]
    first_part <- 0
    second_part <- 0
    
    for (tt in 0:t_max) {
      first_part <- first_part + sum(phi[[tt + 1]] * mu_ext_first[[tt + 1]])
      second_part <- second_part + sum(phi[[tt + 1]] * mu_1_T_0[[tt + 1]])
    }
    
    mu.ext[type] <- first_part - second_part
  }
  
  list(
    Diff_mu = Diff_mu,
    mu.SACE = mu.SACE,
    mu.ext = mu.ext
  )
}
