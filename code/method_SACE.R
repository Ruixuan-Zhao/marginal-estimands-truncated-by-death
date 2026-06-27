## =============================================================================
## --------------------------- Estimate SACE(t) --------------------------------
## =============================================================================

#### Reference: Wang (2017), Biometrika

nLL_beta <- function(beta, gamma, W, sz) {
  
  ebeta <- plogis(W %*% beta)
  egamma <- plogis(W %*% gamma)
  
  loglike <- sum(
    sz[, 1] * log(ebeta) +
      sz[, 2] * log(1 - ebeta) +
      sz[, 3] * log(egamma) +
      sz[, 4] * log(1 - ebeta * egamma)
  )
  
  return(-loglike)
}


nLL_gamma <- function(gamma, beta, W, sz) {
  
  ebeta <- plogis(W %*% beta)
  egamma <- plogis(W %*% gamma)
  
  loglike <- sum(
    sz[, 1] * log(ebeta) +
      sz[, 2] * log(1 - ebeta) +
      sz[, 3] * log(egamma) +
      sz[, 4] * log(1 - ebeta * egamma)
  )
  
  return(-loglike)
}


Alternative <- function(X, A, Z, S, y, beta, gamma) {
  
  lr_fit_z_L0 <- glm(Z ~ A + X, family = binomial())
  p_Z1 <- predict(lr_fit_z_L0, type = "response")
  p_Z0 <- 1 - p_Z1
  
  ### 1. Parameter estimation ##############################################
  
  W <- cbind(rep(1, length(A)), X, A)
  p <- dim(W)[2]
  
  ### Estimating alpha
  W.expit <- (1 - plogis(W %*% gamma)) * Z
  lm.y <- lm(y ~ 1 + X + A + W.expit + Z)
  alpha <- summary(lm.y)$coefficients[, 1]
  
  if (length(alpha) < p + 2) {
    alpha[p + 2] <- 0
  }
  
  ### 2. SACE estimation ###################################################
  
  Z0 <- rep(0, length(A))
  Z1 <- rep(1, length(A))
  
  # 2.1 E[Y(0)|G=LL] equivalent to sace_z0
  
  Coef <- cbind(rep(1, length(A)), X, A, rep(1, length(A)), Z0)
  mu_0_LL_W <- Coef %*% alpha
  mu_0_LL <- mean(mu_0_LL_W * (1 - Z) * S / p_Z0) /
    mean((1 - Z) * S / p_Z0)
  
  # 2.2 E[Y(1)|G=LL] equivalent to sace_z1
  
  Coef <- cbind(rep(1, length(A)), X, A, rep(1, length(A)), Z1)
  mu_1_LL_W <- Coef %*% alpha
  mu_1_LL <- mean(mu_1_LL_W * (1 - Z) * S / p_Z0) /
    mean((1 - Z) * S / p_Z0)
  
  return(mu_1_LL - mu_0_LL)
}


#' Estimate SACE(t)
#'
#' Estimates the survivor average causal effect at each follow-up time using the
#' approach of Wang (2017).
#'
#' @param data A data frame with one row per subject. Required columns are `A`,
#'   the substitution variable; `Z`, the treatment indicator;
#'   baseline covariates named `X1`, `X2`, ...; survival indicators named `S1`,
#'   `S2`, ...; and outcomes named `Y1`, `Y2`, ....
#'
#' @return A named numeric vector containing SACE(t) at each follow-up time.
#'
#' @export
estimate_SACE <- function(data) {
  s_cols <- grep("^S[0-9]+$", names(data), value = TRUE)
  t_max <- length(s_cols)
  x_cols <- grep("^X[0-9]+$", names(data), value = TRUE)
  
  X <- as.matrix(data[x_cols])
  A <- data$A
  Z <- data$Z
  S_all <- as.matrix(data[paste0("S", seq_len(t_max))])
  Y_all <- as.matrix(data[paste0("Y", seq_len(t_max))])
  
  SACE_t <- rep(0, t_max)
  names(SACE_t) <- paste("SACE at time", 1:t_max)
  
  for (t in 1:t_max) {
    S <- S_all[, t]
    
    ### Estimating beta and gamma
    
    s1 <- S == 1
    s0_z1 <- S == 0 & Z == 1
    s1_z0 <- S == 1 & Z == 0
    s0_z0 <- S == 0 & Z == 0
    sz <- cbind(s1, s0_z1, s1_z0, s0_z0)
    
    W <- cbind(rep(1, length(A)), X, A)
    
    thres <- 1e-6
    Diff <- function(x, y) sum((x - y)^2) / sum(x^2 + thres)
    
    beta <- rep(0, ncol(W))
    gamma <- rep(0, ncol(W))
    diff <- thres + 1
    step <- 0
    max.step <- 1000
    
    while (diff > thres & step < max.step) {
      opt1 <- optim(
        beta,
        nLL_beta,
        gamma = gamma,
        W = W,
        sz = sz,
        control = list(maxit = max.step)
      )
      
      diff1 <- Diff(opt1$par, beta)
      beta <- opt1$par
      
      opt2 <- optim(
        gamma,
        nLL_gamma,
        beta = beta,
        W = W,
        sz = sz,
        control = list(maxit = max.step)
      )
      
      diff <- max(diff1, Diff(opt2$par, gamma))
      gamma <- opt2$par
      step <- step + 1
    }
    
    SACE_t[t] <- Alternative(
      X = X,
      A = A,
      Z = Z,
      S = S,
      y = Y_all[, t],
      beta = beta,
      gamma = gamma
    )
  }
  
  SACE_t
}
