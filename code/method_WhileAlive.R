################################################################################
###### Estimation method for while-alive estimands and naive estimands #########
################################################################################

#' Sum components under one weight
#'
#' @keywords internal
#' @noRd
.weighted_lambda_sum <- function(weights, lambda_0_T_0, lambda_1_T_1) {
  lambda <- c("lambda(0)" = 0, "lambda(1)" = 0)

  for (t in seq_along(weights)) {
    lambda["lambda(0)"] <- lambda["lambda(0)"] +
      sum(weights[[t]] * lambda_0_T_0[[t]])
    lambda["lambda(1)"] <- lambda["lambda(1)"] +
      sum(weights[[t]] * lambda_1_T_1[[t]])
  }

  lambda
}

#' Estimate while-alive and naive estimands
#'
#' Estimates while-alive contrasts under the exit-time, average, cumulative,
#' and AUC-based weights, and naive contrasts at each time point.
#'
#' @param data A data frame with one row per subject. Required columns are
#'   `Z`, the treatment indicator; `SURVT`, the observed exit time; `Y0`, the
#'   baseline outcome; and follow-up outcomes `Y1`-`Y[t_max]`, where `t_max`
#'   is the number of follow-up times. All other columns are treated as
#'   baseline covariates.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#'
#' @return A list with:
#'   \describe{
#'     \item{`Diff_lambda`}{While-alive contrasts under the exit-time, average,
#'       cumulative, and AUC-based weights.}
#'     \item{`Diff_lambda.naive`}{naive contrasts at each time point.}
#'   }
#'
#' @export
get_est_WhileAlive_allWeight <- function(data, tau_time) {
  t_max <- length(tau_time) - 1
  outcome_cols <- c("Y0", paste0("Y", seq_len(t_max)))
  baseline_covariate_cols <- setdiff(
    names(data),
    c("Z", "SURVT", outcome_cols)
  )
  baseline_data <- data.frame(
    Z = data$Z,
    data[baseline_covariate_cols],
    SURVT = data$SURVT,
    check.names = FALSE
  )

  survival_prob <- compute_St_Zz(
    data = baseline_data,
    t_max = t_max,
    covariate_cols = baseline_covariate_cols
  )
  weight.sum <- Weights_general(tau_time)
  weights_naive <- Weights_timeSpecified(
    survival_prob$Prob_St_Z0,
    survival_prob$Prob_St_Z1
  )

  tx <- .fit_treatment_prob(
    data = baseline_data,
    covariate_cols = baseline_covariate_cols
  )
  Z <- baseline_data$Z
  p_Z0 <- tx$p_Z0
  p_Z1 <- tx$p_Z1

  Y0 <- data$Y0
  SURVT <- data$SURVT

  lambda_0_T_0 <- .make_components(t_max)
  lambda_1_T_1 <- .make_components(t_max)

  for (t in 0:t_max) {
    lambda_0_T_0[[t + 1]][1] <- mean((Z == 0) * (SURVT == t) * Y0 / p_Z0)
    lambda_1_T_1[[t + 1]][1] <- mean((Z == 1) * (SURVT == t) * Y0 / p_Z1)
  }

  for (r in 1:t_max) {
    y <- data[[paste0("Y", r)]]
    y[is.na(y)] <- 0

    for (t in r:t_max) {
      lambda_0_T_0[[t + 1]][r + 1] <-
        mean((Z == 0) * (SURVT == t) * y / p_Z0)
      lambda_1_T_1[[t + 1]][r + 1] <-
        mean((Z == 1) * (SURVT == t) * y / p_Z1)
    }
  }

  lambda_01.sum <- lapply(
    weight.sum,
    .weighted_lambda_sum,
    lambda_0_T_0,
    lambda_1_T_1
  )

  lambda_01.naive <- vector("list", t_max)
  for (r in 1:t_max) {
    weights_z0 <- weights_naive$weights_timeSpecified_z0[[r]]
    weights_z1 <- weights_naive$weights_timeSpecified_z1[[r]]
    lambda_z0 <- .weighted_lambda_sum(
      weights_z0,
      lambda_0_T_0,
      lambda_1_T_1
    )
    lambda_z1 <- .weighted_lambda_sum(
      weights_z1,
      lambda_0_T_0,
      lambda_1_T_1
    )

    lambda_01.naive[[r]] <- c(
      "lambda(0)" = unname(lambda_z0["lambda(0)"]),
      "lambda(1)" = unname(lambda_z1["lambda(1)"])
    )
  }
  names(lambda_01.naive) <- paste("Naive at time", 1:t_max)

  Diff_lambda <- vapply(
    lambda_01.sum,
    function(lambda) unname(lambda["lambda(1)"] - lambda["lambda(0)"]),
    numeric(1)
  )
  Diff_lambda.naive <- vapply(
    lambda_01.naive,
    function(lambda) unname(lambda["lambda(1)"] - lambda["lambda(0)"]),
    numeric(1)
  )

  list(
    Diff_lambda = Diff_lambda,
    Diff_lambda.naive = Diff_lambda.naive
  )
}
