## ---------------------------
##
## Script name: hyper_beta_omega.R
##
## Purpose of script: elicitates theta on basis of two precision matrices, a threshold c and probability alpha for a tensor spline (adapted from sdPrior)
##
## Author: BLIND
##
## Date Created: 2022-6-5
##
## Email: BLIND
##
## ---------------------------

hyper_beta_omega <- function (K1, K2, alpha = 0.1, ct = 3, R = 10000, myseed = 42) {
  
  omegaseq <- seq(0.05, 0.95, length = 11)
  
  sim_max <- function(scale){
    
    tau2_sample <- rweibull(n = R, shape=0.5, scale=scale)
    omega_sample <- sapply(1:R,  function(x) sample(1:length(omegaseq), size=1))
    
    svd_lookup <- lapply(omegaseq, function(omega) {
      set.seed(myseed)
      # 1) simulate R values of tau2 out of prior
      
      # 2) simulate R beta vector conditional on tau2
      # beta_sample <- sapply(1:R, function(x) rmvnorm(1, sigma=tau2_sample[x] * Kinv, method = "svd"))
      K <- omega * K1  + (1 - omega) * K2
      Kinv <- ginv(K)
      # manual svd
      d <- dim(Kinv)[1]
      s <- svd((Kinv + t(Kinv))/2)
      s$d <- abs(zapsmall(s$d))
      m <- sum(s$d > 0)
      
      # x <- matrix(rnorm(m*R), nrow=m)
      # x <- rbind(x, matrix(0, nrow=d-m, ncol=R))
      # beta_sample <- s$u %*% diag(sqrt(s$d)) %*% x
      # beta_sample
      list(s=s, m=m, d=d)
    })
    
    
    beta_sample <- sapply(omega_sample, function(omega_ind) {
      svd_ind <- svd_lookup[[omega_ind]]
      # x <- matrix(rnorm(m), nrow= svd_lookup[[omega_ind]]$d - m)
      x <- matrix(rnorm(svd_ind$m), nrow=svd_ind$m)
      x <- rbind(x, matrix(0, nrow=svd_ind$d-svd_ind$m, ncol=1))
      svd_ind$s$u %*% diag(sqrt(svd_ind$s$d)) %*% x
    })
    
    beta_sample <- sapply(1:R, function(x) sqrt(tau2_sample[x]) * beta_sample[, x])
    
    
    # max_sample <- apply(abs(beta_sample[-(1:9),]), 1, max)
    max_sample <- apply(abs(beta_sample), 1, max)
    
    max_sample
  }
  
  optfn <- function(scale, ct, alpha){
    quantile(sim_max(scale = scale), probs = alpha) - ct
  }
  
  # 3) determine fraction of examples where max(beta_j) > c for an estimate of the prob
  bopt <- uniroot(f = optfn, interval = c(0.00000001, upper = 1000000), alpha=1-alpha, ct = ct)
  
  bopt
}