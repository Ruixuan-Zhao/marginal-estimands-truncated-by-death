################################################################################
# Shared helper functions for estimand weights and survival probabilities
################################################################################


#' Fit the propensity score model
#'
#' Internal helper used to estimate the treatment probabilities.
#'
#' @param data A data frame containing treatment `Z` and baseline covariates.
#' @param covariate_cols Character vector of baseline covariate column names
#'   used in the model.
#'
#' @return A list with the fitted probabilities `p_Z1 = Pr(Z = 1 | A, X)` and
#'   `p_Z0 = Pr(Z = 0 | A, X)`.
#'
#' @keywords internal
#' @noRd
.fit_treatment_prob <- function(data, covariate_cols) {
  model_data <- data.frame(Z = data$Z, data[covariate_cols], check.names = FALSE)
  fit <- glm(Z ~ ., data = model_data, family = binomial())
  p_Z1 <- predict(fit, type = "response")
  
  list(
    p_Z1 = p_Z1,
    p_Z0 = 1 - p_Z1
  )
}


#' Create empty component list
#'
#' Internal helper that creates the nested list structure used for estimand.
#'
#' @param t_max Maximum follow-up time index.
#'
#' @return A list of length `t_max + 1`. Element `t + 1` is a numeric vector of
#'   zeros with length `t + 1`, corresponding to measurement times
#'   `r = 0, ..., t`.
#'
#' @keywords internal
#' @noRd
.make_components <- function(t_max) {
  components <- vector("list", t_max + 1)

  for (t in 0:t_max) {
    components[[t + 1]] <- rep(0, t + 1)
  }

  components
}


#' Estimate survival probabilities under two treatment arms
#'
#' Computes inverse-probability weighted estimates of `Pr(T(z) >= t)` for
#' `z = 0, 1` and `t = 0, ..., t_max`, where survival through time `t` is encoded
#' as `SURVT >= t`.
#'
#' @param data A data frame containing treatment `Z`, baseline covariates, and 
#'   survival time `SURVT`.
#' @param t_max Maximum follow-up time index.
#' @param covariate_cols Character vector of baseline covariate column names
#'   used in the model.
#'
#' @return A list with two named numeric vectors:
#'   \describe{
#'     \item{`Prob_St_Z0`}{Estimated `Pr(T(0) >= t)` for
#'       `t = 0, ..., t_max`.}
#'     \item{`Prob_St_Z1`}{Estimated `Pr(T(1) >= t)` for
#'       `t = 0, ..., t_max`.}
#'   }
#'
#' @export
compute_St_Zz <- function(data,
                          t_max,
                          covariate_cols) {
  tx <- .fit_treatment_prob(
    data = data,
    covariate_cols = covariate_cols
  )
  Z <- data$Z
  p_Z0 <- tx$p_Z0
  p_Z1 <- tx$p_Z1

  SURVT <- data$SURVT

  Prob_St_Z0 <- rep(1, t_max + 1)
  Prob_St_Z1 <- rep(1, t_max + 1)
  names(Prob_St_Z0) <- 0:t_max
  names(Prob_St_Z1) <- 0:t_max

  for (t in 1:t_max) {
    Prob_St_Z0[t + 1] <- mean((Z == 0) * (SURVT >= t) / p_Z0)
    Prob_St_Z1[t + 1] <- mean((Z == 1) * (SURVT >= t) / p_Z1)
  }

  list(Prob_St_Z0 = Prob_St_Z0, Prob_St_Z1 = Prob_St_Z1)
}


#' Construct general estimand weights
#'
#' Builds the four weighting schemes:
#' exit-time, average, cumulative, and AUC-based weights. For each time point
#' `t`, the corresponding vector contains weights for
#' `r = 0, ..., t`.
#'
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#'
#' @return A named list with elements `Exit time`, `Average`, `Cumulative`, and
#'   `AUC-based`.
#'
#' @export
Weights_general <- function(tau_time) {
  t_max <- length(tau_time) - 1

  weight_exit <- vector("list", t_max + 1)
  weight_ave <- vector("list", t_max + 1)
  weight_cum <- vector("list", t_max + 1)
  weight_auc <- vector("list", t_max + 1)

  for (t in 0:t_max) {
    n_time_points <- t + 1

    weight_exit[[t + 1]] <- rep(0, n_time_points)
    weight_exit[[t + 1]][n_time_points] <- 1

    weight_ave[[t + 1]] <- rep(1 / n_time_points, n_time_points)
    weight_cum[[t + 1]] <- rep(1, n_time_points)

    weight_auc[[t + 1]] <- numeric(n_time_points)
    if (t > 0) {
      for (r in 0:t) {
        weight_auc[[t + 1]][r + 1] <- if (r == 0) {
          0.5 * tau_time[2]
        } else if (r == t) {
          0.5 * (tau_time[t + 1] - tau_time[t])
        } else {
          0.5 * (tau_time[r + 2] - tau_time[r])
        }
      }
    }
  }

  weight.sum <- list(weight_exit, weight_ave, weight_cum, weight_auc)
  names(weight.sum) <- c("Exit time", "Average", "Cumulative", "AUC-based")

  weight.sum
}


#' Construct time-specific weights
#'
#' @param Prob_St_Z0 Numeric vector of `Pr(T(0) >= t)` for `t = 0, ..., t_max`.
#' @param Prob_St_Z1 Numeric vector of `Pr(T(1) >= t)` for `t = 0, ..., t_max`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`weights_timeSpecified_z0`}{A list of time-specific weight objects
#'       using survival probabilities for `Z = 0`.}
#'     \item{`weights_timeSpecified_z1`}{A list of time-specific weight objects
#'       using survival probabilities for `Z = 1`.}
#'   }
#'
#' @keywords internal
#' @noRd
Weights_timeSpecified <- function(Prob_St_Z0, Prob_St_Z1) {
  t_max <- length(Prob_St_Z0) - 1

  weights_timeSpecified_z0 <- vector("list", t_max)
  weights_timeSpecified_z1 <- vector("list", t_max)

  for (r in 1:t_max) {
    weight_z0 <- .make_components(t_max)
    weight_z1 <- .make_components(t_max)

    for (t in r:t_max) {
      weight_z0[[t + 1]][r + 1] <- 1 / Prob_St_Z0[r + 1]
      weight_z1[[t + 1]][r + 1] <- 1 / Prob_St_Z1[r + 1]
    }

    weights_timeSpecified_z0[[r]] <- weight_z0
    weights_timeSpecified_z1[[r]] <- weight_z1
  }

  list(
    weights_timeSpecified_z0 = weights_timeSpecified_z0,
    weights_timeSpecified_z1 = weights_timeSpecified_z1
  )
}


#' Clamp probabilities away from 0 and 1
#'
#' @keywords internal
#' @noRd
.clamp_prob <- function(p, eps = 1e-8) {
  pmin(pmax(p, eps), 1 - eps)
}
