################################################################################
###### Generate simulated data for marginal separable effect estimands #########
################################################################################

#' Generate simulated data for marginal separable effect estimands
#'
#' @param n.sample Number of samples.
#' @param seed Optional random seed.
#'
#' @return a data frame containing the observed simulated data. 
#'   Columns include baseline covariates `A` and `X1`-`X3`, treatment `Z`, 
#'   survival indicators `S1`-`S3`, 
#'   observed survival time `SURVT`, time-varying covariates `L1`-`L3`, 
#'   baseline outcome `Y0`, and follow-up outcomes`Y1`-`Y3`.
#'
#' @export
simulate_SE <- function(n.sample, seed) {
  if (!is.null(seed)) set.seed(seed)

  # Baseline covariates L0 = (X, A).
  X1 <- sample(c(-1, 1), size = n.sample, replace = TRUE)
  X2 <- runif(n.sample, min = -1, max = 1)
  X3 <- runif(n.sample, min = -1, max = 1)

  # Substitution variable.
  p_A <- plogis(0.0 + 0.2 * X1 + 0.1 * X2 - 0.1 * X3)
  A <- rbinom(n.sample, size = 1, prob = p_A)
  L0 <- cbind(X1, X2, X3, A)

  # Exposure components. Set Zy = Zs = Z.
  p_Zy <- plogis(0.0 + 0.1 * A + 0.2 * X1 - 0.1 * X2 + 0.1 * X3)
  Zy <- rbinom(n.sample, size = 1, prob = p_Zy)
  Zs <- Zy
  Z <- Zy
  
  unit1 <- c(1 / 2, 1 / 2, 1 / 2, 1)
  survival_shift <- -1.1
  survival_Zs_coef <- 0.5

  # Time 1.
  p_L1 <- plogis(L0 %*% unit1 + Zs)
  L1 <- rbinom(n.sample, size = 1, prob = p_L1)
  p_S1 <- plogis(survival_shift + L0 %*% unit1 + L1 + survival_Zs_coef * Zs)
  S1 <- rbinom(n.sample, size = 1, prob = p_S1)
  Y1 <- rep(NA, n.sample)
  mean_Y1 <- 1 / 2 + L0 %*% unit1 + L1 + Zy
  Y1[S1 == 1] <- rnorm(sum(S1), mean = mean_Y1[S1 == 1], sd = 0.5)

  # Time 2.
  L2 <- L1
  p_L2 <- plogis(L0 %*% unit1 + L1 + Zs)
  L2[S1 == 1] <- rbinom(sum(S1), size = 1, prob = p_L2[S1 == 1])
  S2 <- rep(0, n.sample)
  p_S2 <- plogis(survival_shift + L0 %*% unit1 + L2 + survival_Zs_coef * Zs)
  S2[S1 == 1] <- rbinom(sum(S1), size = 1, prob = p_S2[S1 == 1])
  Y2 <- rep(NA, n.sample)
  mean_Y2 <- 1 / 2 + L0 %*% unit1 + L2 + Zy + Y1
  Y2[S2 == 1] <- rnorm(sum(S2), mean = mean_Y2[S2 == 1], sd = 0.5)

  # Time 3.
  L3 <- L2
  p_L3 <- plogis(L0 %*% unit1 + L2 + Zs)
  L3[S2 == 1] <- rbinom(sum(S2), size = 1, prob = p_L3[S2 == 1])
  S3 <- rep(0, n.sample)
  p_S3 <- plogis(survival_shift + L0 %*% unit1 + L3 + survival_Zs_coef * Zs)
  S3[S2 == 1] <- rbinom(sum(S2), size = 1, prob = p_S3[S2 == 1])
  Y3 <- rep(NA, n.sample)
  mean_Y3 <- 1 / 2 + L0 %*% unit1 + L3 + Zy + Y2
  Y3[S3 == 1] <- rnorm(sum(S3), mean = mean_Y3[S3 == 1], sd = 0.5)

  # Survival time and baseline outcome.
  SURVT <- apply(cbind(S1, S2, S3), 1, sum)
  Y0 <- rep(0, n.sample)

  data <- data.frame(
    A = A,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    Z = Z,
    S1 = S1,
    S2 = S2,
    S3 = S3,
    SURVT = SURVT,
    L1 = L1,
    L2 = L2,
    L3 = L3,
    Y0 = Y0,
    Y1 = Y1,
    Y2 = Y2,
    Y3 = Y3
  )

  data
}
