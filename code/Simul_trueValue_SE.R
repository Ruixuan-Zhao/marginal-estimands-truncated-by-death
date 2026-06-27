################################################################################
######## Simulation: True value for marginal separable effect estimands ########
################################################################################


Prop_L_func <- function(Lr, Lr_1, L0, Z) {
  unit1 <- c(1 / 2, 1 / 2, 1 / 2, 1)
  eta <- if (is.na(Lr_1)) {
    L0 %*% unit1 + Z
  } else {
    L0 %*% unit1 + Lr_1 + Z
  }
  
  p <- plogis(eta)
  if (Lr == 1) p else 1 - p
}



Prop_S_func <- function(Sr, L0, Lr, Z) {
  unit1 <- c(1 / 2, 1 / 2, 1 / 2, 1)
  survival_shift <- -1.1
  survival_Zs_coef <- 0.5
  
  p <- plogis(survival_shift + L0 %*% unit1 + Lr + survival_Zs_coef * Z)
  if (Sr == 1) p else 1 - p
}


#' Compute truth values for marginal separable effect estimands in simulation study
#'
#' @param n.sample Number of samples.
#' @param seed Optional random seed.
#'
#' @return A list with:
#'   \describe{
#'     \item{`trueValue_Prob_T`}{A matrix of \eqn{\Pr(T = t \mid Z = z_s)} for
#'       \eqn{t = 0, 1, 2, 3} and \eqn{z_s = 0, 1}.}
#'     \item{`trueValue_SE`}{A matrix of truth values for the marginal
#'       separable effect contrasts \eqn{\Gamma(1, z_S) - \Gamma(0, z_S)}.}
#'   }
#'
#' @export
compute_truth_values_SE <- function(n.sample, seed) {
  sim_data <- simulate_SE(n.sample = n.sample, seed = seed)
  L0 <- as.matrix(sim_data[c("X1", "X2", "X3", "A")])

  trueValue_Prob_T <- matrix(0, 4, 2)
  colnames(trueValue_Prob_T) <- c("z_s=0", "z_s=1")
  rownames(trueValue_Prob_T) <- c("t=0", "t=1", "t=2", "t=3")

  # Pr(T = 0 | Z = z_s, L0).
  Prob_T0_z0 <- Prop_S_func(Sr = 0, L0 = L0, Lr = 1, Z = 0) *
    Prop_L_func(Lr = 1, Lr_1 = NA, L0 = L0, Z = 0) +
    Prop_S_func(Sr = 0, L0 = L0, Lr = 0, Z = 0) *
    Prop_L_func(Lr = 0, Lr_1 = NA, L0 = L0, Z = 0)
  Prob_T0_z1 <- Prop_S_func(Sr = 0, L0 = L0, Lr = 1, Z = 1) *
    Prop_L_func(Lr = 1, Lr_1 = NA, L0 = L0, Z = 1) +
    Prop_S_func(Sr = 0, L0 = L0, Lr = 0, Z = 1) *
    Prop_L_func(Lr = 0, Lr_1 = NA, L0 = L0, Z = 1)

  trueValue_Prob_T[1, 1] <- mean(Prob_T0_z0)
  trueValue_Prob_T[1, 2] <- mean(Prob_T0_z1)

  # Pr(T = 1 | Z = z_s, L0).
  L2_tilde <- expand.grid(c(0, 1), c(0, 1))
  Prob_T1_z0 <- 0
  Prob_T1_z1 <- 0

  for (i in 1:nrow(L2_tilde)) {
    Prob_T1_z0 <- Prob_T1_z0 +
      Prop_S_func(Sr = 0, L0 = L0, Lr = L2_tilde[i, 2], Z = 0) *
      Prop_L_func(Lr = L2_tilde[i, 2], Lr_1 = L2_tilde[i, 1], L0 = L0, Z = 0) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L2_tilde[i, 1], Z = 0) *
      Prop_L_func(Lr = L2_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 0)

    Prob_T1_z1 <- Prob_T1_z1 +
      Prop_S_func(Sr = 0, L0 = L0, Lr = L2_tilde[i, 2], Z = 1) *
      Prop_L_func(Lr = L2_tilde[i, 2], Lr_1 = L2_tilde[i, 1], L0 = L0, Z = 1) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L2_tilde[i, 1], Z = 1) *
      Prop_L_func(Lr = L2_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 1)
  }

  trueValue_Prob_T[2, 1] <- mean(Prob_T1_z0)
  trueValue_Prob_T[2, 2] <- mean(Prob_T1_z1)

  # Pr(T = 2 | Z = z_s, L0).
  L3_tilde <- expand.grid(c(0, 1), c(0, 1), c(0, 1))
  Prob_T2_z0 <- 0
  Prob_T2_z1 <- 0

  for (i in 1:nrow(L3_tilde)) {
    Prob_T2_z0 <- Prob_T2_z0 +
      Prop_S_func(Sr = 0, L0 = L0, Lr = L3_tilde[i, 3], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 3], Lr_1 = L3_tilde[i, 2], L0 = L0, Z = 0) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 2], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 2], Lr_1 = L3_tilde[i, 1], L0 = L0, Z = 0) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 1], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 0)

    Prob_T2_z1 <- Prob_T2_z1 +
      Prop_S_func(Sr = 0, L0 = L0, Lr = L3_tilde[i, 3], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 3], Lr_1 = L3_tilde[i, 2], L0 = L0, Z = 1) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 2], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 2], Lr_1 = L3_tilde[i, 1], L0 = L0, Z = 1) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 1], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 1)
  }

  trueValue_Prob_T[3, 1] <- mean(Prob_T2_z0)
  trueValue_Prob_T[3, 2] <- mean(Prob_T2_z1)

  # Pr(T = 3 | Z = z_s, L0).
  Prob_T3_z0 <- 0
  Prob_T3_z1 <- 0

  for (i in 1:nrow(L3_tilde)) {
    Prob_T3_z0 <- Prob_T3_z0 +
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 3], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 3], Lr_1 = L3_tilde[i, 2], L0 = L0, Z = 0) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 2], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 2], Lr_1 = L3_tilde[i, 1], L0 = L0, Z = 0) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 1], Z = 0) *
      Prop_L_func(Lr = L3_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 0)

    Prob_T3_z1 <- Prob_T3_z1 +
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 3], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 3], Lr_1 = L3_tilde[i, 2], L0 = L0, Z = 1) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 2], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 2], Lr_1 = L3_tilde[i, 1], L0 = L0, Z = 1) *
      Prop_S_func(Sr = 1, L0 = L0, Lr = L3_tilde[i, 1], Z = 1) *
      Prop_L_func(Lr = L3_tilde[i, 1], Lr_1 = NA, L0 = L0, Z = 1)
  }

  trueValue_Prob_T[4, 1] <- mean(Prob_T3_z0)
  trueValue_Prob_T[4, 2] <- mean(Prob_T3_z1)

  trueValue_SE <- matrix(0, 4, 2)
  colnames(trueValue_SE) <- c(
    "Gamma(1,0)-Gamma(0,0)(z_s=0)",
    "Gamma(1,1)-Gamma(0,1)(z_s=1)"
  )
  rownames(trueValue_SE) <- c("Exit time", "Average", "Cumulative", "AUC-based")

  trueValue_SE[1, 1] <- sum(trueValue_Prob_T[2:4, 1] * c(1, 2, 3))
  trueValue_SE[1, 2] <- sum(trueValue_Prob_T[2:4, 2] * c(1, 2, 3))
  trueValue_SE[2, 1] <- sum(trueValue_Prob_T[2:4, 1] * c(0.5, 1, 1.5))
  trueValue_SE[2, 2] <- sum(trueValue_Prob_T[2:4, 2] * c(0.5, 1, 1.5))
  trueValue_SE[3, 1] <- sum(trueValue_Prob_T[2:4, 1] * c(1, 3, 6))
  trueValue_SE[3, 2] <- sum(trueValue_Prob_T[2:4, 2] * c(1, 3, 6))
  trueValue_SE[4, 1] <- sum(trueValue_Prob_T[2:4, 1] * c(1 / 8, 1 / 2, 7 / 4))
  trueValue_SE[4, 2] <- sum(trueValue_Prob_T[2:4, 2] * c(1 / 8, 1 / 2, 7 / 4))

  list(
    trueValue_Prob_T = trueValue_Prob_T,
    trueValue_SE = trueValue_SE
  )
}
