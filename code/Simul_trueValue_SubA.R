################################################################################
########################### Simulation: True value for #########################
###### While guaranteed-survival and while extended-survival estimands #########
################################################################################

# ============================================================
# Model-based probabilities for Pr(T(0) = t)
# ============================================================

#' Compute `Pr(T(0) = t)` for `t = 0, 1, 2, 3`
#'
#' @param potential_data Model components returned by `simulate_SubA()`.
#'
#' @return A named numeric vector with probabilities `prob_T0 = Pr(T(0) = t)` 
#'   for `t = 0, 1, 2, 3`.
#'
#' @keywords internal
#' @noRd
calc_prob_T0_from_model <- function(potential_data) {
  
  p1rho1 <- potential_data$p1 * potential_data$rho1
  p2rho2 <- potential_data$p2 * potential_data$rho2
  p3rho3 <- potential_data$p3 * potential_data$rho3
  
  prob_T0_0 <- mean(1 - p1rho1)
  
  prob_T0_1 <- mean(
    p1rho1 * (1 - p2rho2)
  )
  
  prob_T0_2 <- mean(
    p1rho1 * p2rho2 * (1 - p3rho3)
  )
  
  prob_T0_3 <- mean(
    p1rho1 * p2rho2 * p3rho3
  )
  
  prob_T0 <- c(prob_T0_0, prob_T0_1, prob_T0_2, prob_T0_3)
  names(prob_T0) <- paste0("T0_", 0:3)
  
  return(prob_T0)
}


# ============================================================
# Truth for mu(1) - mu(0)
#
# Under the current DGP and Y^0 = 0:
#
# mu(1) - mu(0)
# = sum_{t=1}^3 sum_{r=1}^t w_t^r * (eta_Z + eta_G) * r * Pr{T(0)=t}
#
# Here Pr{T(0)=t} is computed from the model probabilities
# p_i(L^0) and rho_i(L^0), not from realized T0.
# ============================================================

#' Compute truth value for \eqn{\mu(1) - \mu(0)} under one weighting 
#' scheme
#'
#' @param potential_data Model components returned by `simulate_SubA()`.
#' @param weight.list A nested weight list for one weighting scheme.
#'
#' @return A named numeric scalar for \eqn{\mu(1) - \mu(0)}.
#'
#' @keywords internal
#' @noRd
calc_truth_mu_diff_one_weight <- function(potential_data, weight.list) {
  
  t_max <- length(weight.list) - 1
  
  prob_T0 <- calc_prob_T0_from_model(potential_data)
  
  eta_G <- potential_data$eta_G[1]
  eta_Z <- potential_data$eta_Z[1]
  
  mu_diff <- 0
  
  for (t in 1:t_max) {
    
    prob_T0_t <- prob_T0[t + 1]
    
    for (r in 1:t) {
      
      w_tr <- weight.list[[t + 1]][r + 1]
      
      mu_diff <- mu_diff + w_tr * (eta_Z + eta_G) * r * prob_T0_t
    }
  }
  
  return(c(mu_diff = as.numeric(mu_diff)))
}


# ============================================================
# Model-based components for mu_ext
#
# mu_ext
# = sum_{t=1}^3 sum_{r=1}^t phi_t^r E[
#   {a_t(L0)-b_t(L0)}
#   {r B0(A,X) + sum_{k=1}^r [1 - P_B(k|L0)]}
#   + eta_G * r {a_t(L0) P_K(r|L0) - b_t(L0)}
# ]
#
# where
# a_t(L0) = Pr{T(1)=t | L0}
# b_t(L0) = Pr{T(0)=t | L0}
# P_B(k|L0) = Pr(T_B > k | L0)
# P_K(r|L0) = Pr(T_K > r | L0)
# ============================================================

#' Compute one weighted truth value for the while extended-survival 
#' estimand \eqn{\mu^{\text{ext}}}
#'
#' @param potential_data Model components returned by `simulate_SubA()`.
#' @param phi.list A nested weight list for one weighting scheme.
#'
#' @return A named numeric scalar for \eqn{\mu^{\text{ext}}}.
#'
#' @keywords internal
#' @noRd
calc_truth_mu_ext_one_weight <- function(potential_data, phi.list) {
  
  t_max <- length(phi.list) - 1
  n <- nrow(potential_data)
  
  eta_G <- potential_data$eta_G[1]
  eta_Z <- potential_data$eta_Z[1]
  
  # P_B(k | L0) = Pr(T_B > k | L0)
  PB <- list()
  PB[[1]] <- rep(1, n)  # k = 0
  PB[[2]] <- potential_data$s1
  PB[[3]] <- potential_data$s1 * potential_data$s2
  PB[[4]] <- potential_data$s1 * potential_data$s2 * potential_data$s3
  
  # P_K(r | L0) = Pr(T_K > r | L0)
  PK <- list()
  PK[[1]] <- rep(1, n)  # r = 0
  PK[[2]] <- potential_data$rho1
  PK[[3]] <- potential_data$rho1 * potential_data$rho2
  PK[[4]] <- potential_data$rho1 * potential_data$rho2 * potential_data$rho3
  
  # a_t(L0) = Pr(T(1)=t | L0)
  a <- list()
  a[[1]] <- NULL
  a[[2]] <- potential_data$p1 * (1 - potential_data$p2)
  a[[3]] <- potential_data$p1 * potential_data$p2 * (1 - potential_data$p3)
  a[[4]] <- potential_data$p1 * potential_data$p2 * potential_data$p3
  
  # b_t(L0) = Pr(T(0)=t | L0)
  p1rho1 <- potential_data$p1 * potential_data$rho1
  p2rho2 <- potential_data$p2 * potential_data$rho2
  p3rho3 <- potential_data$p3 * potential_data$rho3
  
  b <- list()
  b[[1]] <- NULL
  b[[2]] <- p1rho1 * (1 - p2rho2)
  b[[3]] <- p1rho1 * p2rho2 * (1 - p3rho3)
  b[[4]] <- p1rho1 * p2rho2 * p3rho3
  
  mu_ext <- 0
  
  for (t in 1:t_max) {
    
    a_t <- a[[t + 1]]
    b_t <- b[[t + 1]]
    
    for (r in 1:t) {
      
      phi_tr <- phi.list[[t + 1]][r + 1]
      
      # sum_{k=1}^r {1 - P_B(k | L0)}
      sum_progression <- rep(0, n)
      for (k in 1:r) {
        sum_progression <- sum_progression + (1 - PB[[k + 1]])
      }
      
      base_part <- r * (potential_data$B0 + eta_Z) + sum_progression
      
      PK_r <- PK[[r + 1]]
      
      term_ext <- mean(
        (a_t - b_t) * base_part +
          eta_G * r * (a_t * PK_r - b_t)
      )
      
      mu_ext <- mu_ext + phi_tr * term_ext
    }
  }
  
  return(c(mu_ext = as.numeric(mu_ext)))
}


# ============================================================
# Compute all truth values
# ============================================================

#' Compute truth values for the While guaranteed-survival and while
#' extended-survival estimands in simulation study
#'
#' Computes the simulation truth values for the while guaranteed-survival and
#' while extended-survival estimands using model components returned by
#' `simulate_SubA(..., keep_potential = TRUE)`.
#'
#' @param potential_data A data frame containing the model components needed for
#'   truth calculation: `B0`, `p1`-`p3`, `rho1`-`rho3`, `s1`-`s3`, `eta_G`,
#'   and `eta_Z`.
#' @param tau_time Numeric vector of follow-up times, starting at baseline.
#'
#' @return A list with two data frames:
#'   \describe{
#'     \item{`mu_diff`}{Truth values for \eqn{\mu(1) - \mu(0)} under
#'       the exit-time, average, cumulative, and AUC-based weights.}
#'     \item{`mu_ext`}{Truth values for \eqn{\mu^{\text{ext}}} under
#'       the cumulative and AUC-based weights.}
#'   }
#'
#' @export
compute_truth_values <- function(potential_data,
                                 tau_time = c(0, 1/4, 1/2, 1)) {
  
  weight.sum <- Weights_general(tau_time)
  
  # ------------------------------------------------------------
  # mu(1) - mu(0)
  # Use all four weight choices:
  # Exit time, Average, Cumulative, AUC-based
  # ------------------------------------------------------------
  
  mu_diff_results <- do.call(
    rbind,
    lapply(seq_along(weight.sum), function(j) {
      
      val <- calc_truth_mu_diff_one_weight(
        potential_data = potential_data,
        weight.list = weight.sum[[j]]
      )
      
      data.frame(
        weight_type = names(weight.sum)[j],
        mu_diff = val["mu_diff"],
        row.names = NULL
      )
    })
  )
  
  # ------------------------------------------------------------
  # mu_ext
  # Use only Cumulative and AUC-based weights for phi
  # ------------------------------------------------------------
  
  phi_indices <- c(3, 4)
  
  mu_ext_results <- do.call(
    rbind,
    lapply(phi_indices, function(j) {
      
      val <- calc_truth_mu_ext_one_weight(
        potential_data = potential_data,
        phi.list = weight.sum[[j]]
      )
      
      data.frame(
        phi_type = names(weight.sum)[j],
        mu_ext = val["mu_ext"],
        row.names = NULL
      )
    })
  )
  
  return(list(
    mu_diff = mu_diff_results,
    mu_ext = mu_ext_results
  ))
}
