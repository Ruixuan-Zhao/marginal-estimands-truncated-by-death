################################################################################
######################## Generate simulated data ###############################
###### While guaranteed-survival and while extended-survival estimands #########
################################################################################

#' Generate simulated data for the While guaranteed-survival and while
#' extended-survival estimands in simulation study
#'
#' @param n.sample Number of samples.
#' @param seed Optional random seed. 
#' @param keep_potential Logical. If `TRUE`, return the latent
#'   data components used to compute simulation truth values.
#'
#' @return A list with element `data`, a data frame containing the observed
#'   simulated data. If `keep_potential = TRUE`, the list also contains
#'   `potential_data`, a data frame with model components used by
#'   `compute_truth_values()`.
#'
#' @export
simulate_SubA <- function(n.sample,
                          seed = NULL,
                          keep_potential = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # ============================================================
  # 1. Baseline covariates
  # ============================================================
  
  X1 <- sample(c(-1, 1), size = n.sample, replace = TRUE)
  X2 <- runif(n.sample, min = -1, max = 1)
  X3 <- runif(n.sample, min = -1, max = 1)
  
  pA <- plogis(0.0 + 0.2 * X1 + 0.1 * X2 - 0.1 * X3)
  A <- rbinom(n.sample, size = 1, prob = pA)
  
  X <- cbind(X1, X2, X3)
  colnames(X) <- c("X1", "X2", "X3")
  
  # ============================================================
  # 1.5 Linear predictors
  # ============================================================
  # H determines treatment-resistant lethal survival R_C.
  # R determines treatment-killable harmful survival R_K.
  # D determines disease progression R_B.
  #
  # W order for beta/gamma:
  # W = (1, X1, X2, X3, A)
  # ============================================================
  
  H <- 0.2 * A + 0.3 * X1 - 0.2 * X2 + 0.1 * X3
  R <- 0.1 * A - 0.2 * X1 + 0.1 * X2 + 0.2 * X3
  D <- 0.2 * A - 0.1 * X1 + 0.2 * X2 - 0.1 * X3
  
  # ============================================================
  # 2. Treatment assignment
  # ============================================================
  e <- plogis(0.0 + 0.1 * A + 0.2 * X1 - 0.1 * X2 + 0.1 * X3)
  Z <- rbinom(n.sample, size = 1, prob = e)
  
  # ============================================================
  # 3. Latent mechanism uniforms
  # ============================================================
  
  UK1 <- runif(n.sample)
  UK2 <- runif(n.sample)
  UK3 <- runif(n.sample)
  
  UB1 <- runif(n.sample)
  UB2 <- runif(n.sample)
  UB3 <- runif(n.sample)
  
  UC1 <- runif(n.sample)
  UC2 <- runif(n.sample)
  UC3 <- runif(n.sample)
  
  # ============================================================
  # 4. Treatment-resistant lethal mechanism T_C
  # R_C^r = I(T_C > r)
  #
  # True beta values, W = (1, X1, X2, X3, A):
  # beta1 = c(2.2, 0.3, -0.2, 0.1, 0.2)
  # beta2 = c(2.1, 0.3, -0.2, 0.1, 0.2)
  # beta3 = c(2.0, 0.3, -0.2, 0.1, 0.2)
  # ============================================================
  
  p1 <- plogis(2.2 + H)
  p2 <- plogis(2.1 + H)
  p3 <- plogis(2.0 + H)
  
  RC0 <- rep(1L, n.sample)
  RC1 <- RC0 * as.integer(UC1 <= p1)
  RC2 <- RC1 * as.integer(UC2 <= p2)
  RC3 <- RC2 * as.integer(UC3 <= p3)
  
  TC <- rep(Inf, n.sample)
  TC[RC1 == 0] <- 1
  TC[RC1 == 1 & RC2 == 0] <- 2
  TC[RC2 == 1 & RC3 == 0] <- 3
  
  # ============================================================
  # 5. Treatment-killable harmful mechanism T_K
  # R_K^r = I(T_K > r)
  #
  # True gamma values, W = (1, X1, X2, X3, A):
  # gamma1 = c(1.4, -0.2, 0.1, 0.2, 0.1)
  # gamma2 = c(1.4, -0.2, 0.1, 0.2, 0.1)
  # gamma3 = c(1.4, -0.2, 0.1, 0.2, 0.1)
  #
  # This gamma controls the survival ratio:
  # Pr(S^r(0)=1 | S^{r-1}(0)=1, L0)
  # ------------------------------------------------------------
  # Pr(S^r(1)=1 | S^{r-1}(1)=1, L0)
  # ============================================================
  
  rho1 <- plogis(1.4 + R)
  rho2 <- plogis(1.4 + R)
  rho3 <- plogis(1.4 + R)
  
  RK0 <- rep(1L, n.sample)
  RK1 <- RK0 * as.integer(UK1 <= rho1)
  RK2 <- RK1 * as.integer(UK2 <= rho2)
  RK3 <- RK2 * as.integer(UK3 <= rho3)
  
  TK <- rep(Inf, n.sample)
  TK[RK1 == 0] <- 1
  TK[RK1 == 1 & RK2 == 0] <- 2
  TK[RK2 == 1 & RK3 == 0] <- 3
  
  # ============================================================
  # 6. Treatment-resistant disease progression mechanism T_B
  # R_B^r = I(T_B > r)
  # ============================================================
  
  s1 <- plogis(0.5 - D)
  s2 <- plogis(0.4 - D)
  s3 <- plogis(0.3 - D)
  
  RB0 <- rep(1L, n.sample)
  RB1 <- RB0 * as.integer(UB1 <= s1)
  RB2 <- RB1 * as.integer(UB2 <= s2)
  RB3 <- RB2 * as.integer(UB3 <= s3)
  
  TB <- rep(Inf, n.sample)
  TB[RB1 == 0] <- 1
  TB[RB1 == 1 & RB2 == 0] <- 2
  TB[RB2 == 1 & RB3 == 0] <- 3
  
  # ============================================================
  # 7. Potential survival
  # S^r(0) = I(T_K > r, T_C > r)
  # S^r(1) = I(T_C > r)
  # ============================================================
  
  S1_0 <- RK1 * RC1
  S2_0 <- RK2 * RC2
  S3_0 <- RK3 * RC3
  
  S1_1 <- RC1
  S2_1 <- RC2
  S3_1 <- RC3
  
  T0 <- S1_0 + S2_0 + S3_0
  T1 <- S1_1 + S2_1 + S3_1
  
  # ============================================================
  # 8. Potential disease progression status
  # L^r(0) = I(T_K <= r or T_B <= r or T_C <= r)
  # L^r(1) = I(T_B <= r or T_C <= r)
  # ============================================================
  
  L1_0 <- as.integer(RK1 == 0 | RB1 == 0 | RC1 == 0)
  L2_0 <- as.integer(RK2 == 0 | RB2 == 0 | RC2 == 0)
  L3_0 <- as.integer(RK3 == 0 | RB3 == 0 | RC3 == 0)
  
  L1_1 <- as.integer(RB1 == 0 | RC1 == 0)
  L2_1 <- as.integer(RB2 == 0 | RC2 == 0)
  L3_1 <- as.integer(RB3 == 0 | RC3 == 0)
  
  # ============================================================
  # 9. Principal strata
  # G^r = (S^r(1), S^r(0))
  # ============================================================
  
  make_G <- function(S_1, S_0) {
    out <- rep(NA_character_, length(S_1))
    out[S_1 == 1 & S_0 == 1] <- "LL"
    out[S_1 == 1 & S_0 == 0] <- "LD"
    out[S_1 == 0 & S_0 == 0] <- "DD"
    out[S_1 == 0 & S_0 == 1] <- "DL"
    out
  }
  
  G1 <- make_G(S1_1, S1_0)
  G2 <- make_G(S2_1, S2_0)
  G3 <- make_G(S3_1, S3_0)
  
  G1_num <- as.integer(G1 == "LL")
  G2_num <- as.integer(G2 == "LL")
  G3_num <- as.integer(G3 == "LL")
  
  # ============================================================
  # 10. Potential outcomes
  # ============================================================
  
  Bbase <- 0.5 + 0.2 * A + 0.3 * X1 - 0.2 * X2 + 0.2 * X3
  
  eps1_0 <- rnorm(n.sample, mean = 0, sd = 0.5)
  eps2_0 <- rnorm(n.sample, mean = 0, sd = 0.5)
  eps3_0 <- rnorm(n.sample, mean = 0, sd = 0.5)
  
  eps1_1 <- rnorm(n.sample, mean = 0, sd = 0.5)
  eps2_1 <- rnorm(n.sample, mean = 0, sd = 0.5)
  eps3_1 <- rnorm(n.sample, mean = 0, sd = 0.5)
  
  Y0_0 <- rep(0, n.sample)
  Y0_1 <- rep(0, n.sample)
  eta_G <- 0.1
  eta_Z <- 0.5
  
  # Potential outcomes under treatment
  Y1_1 <- ifelse(
    S1_1 == 1,
    Bbase + L1_1 + eta_Z + eta_G * G1_num + eps1_1,
    NA_real_
  )
  
  Y2_1 <- ifelse(
    S2_1 == 1,
    Bbase + L2_1 + eta_Z + 2 * eta_G * G2_num + Y1_1 - eta_G * G1_num + eps2_1,
    NA_real_
  )
  
  Y3_1 <- ifelse(
    S3_1 == 1,
    Bbase + L3_1 + eta_Z + 3 * eta_G * G3_num + Y2_1 - 2 * eta_G * G2_num + eps3_1,
    NA_real_
  )
  
  # Potential outcome under control
  Y1_0 <- ifelse(
    S1_0 == 1,
    Bbase + L1_0 + eps1_0,
    NA_real_
  )
  
  Y2_0 <- ifelse(
    S2_0 == 1,
    Bbase + L2_0 + Y1_0 + eps2_0,
    NA_real_
  )
  
  Y3_0 <- ifelse(
    S3_0 == 1,
    Bbase + L3_0 + Y2_0 + eps3_0,
    NA_real_
  )
  
  # ============================================================
  # 11. Observed data
  # ============================================================
  
  S1 <- ifelse(Z == 1, S1_1, S1_0)
  S2 <- ifelse(Z == 1, S2_1, S2_0)
  S3 <- ifelse(Z == 1, S3_1, S3_0)
  
  SURVT <- S1 + S2 + S3
  
  L1 <- ifelse(S1 == 1, ifelse(Z == 1, L1_1, L1_0), NA_integer_)
  L2 <- ifelse(S2 == 1, ifelse(Z == 1, L2_1, L2_0), NA_integer_)
  L3 <- ifelse(S3 == 1, ifelse(Z == 1, L3_1, L3_0), NA_integer_)
  
  Y1 <- ifelse(S1 == 1, ifelse(Z == 1, Y1_1, Y1_0), NA_real_)
  Y2 <- ifelse(S2 == 1, ifelse(Z == 1, Y2_1, Y2_0), NA_real_)
  Y3 <- ifelse(S3 == 1, ifelse(Z == 1, Y3_1, Y3_0), NA_real_)
  
  # ============================================================
  # 12. Observed data
  # ============================================================
  
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
    
    Y0 = rep(0, n.sample),
    Y1 = Y1,
    Y2 = Y2,
    Y3 = Y3
  )
  
  out <- list(data = data)
  
  if (keep_potential) {
    potential_data <- data.frame(
      B0 = Bbase,
      
      p1 = p1,
      p2 = p2,
      p3 = p3,
      
      rho1 = rho1,
      rho2 = rho2,
      rho3 = rho3,
      
      s1 = s1,
      s2 = s2,
      s3 = s3,
      
      eta_G = eta_G,
      eta_Z = eta_Z
    )
    
    out$potential_data <- potential_data
  }
  
  return(out)
}
