##
## Script name: 04_fram_te_sd.R
##
## Purpose of script: estimates model bctm_tensor_sd for Framingham
##
## Author: BLIND
##
## Date Created: 2022-10-7
##
## Email: BLIND
##
## ---------------------------

library(dplyr)
# library("RhpcBLASctl")
# omp_set_num_threads(1)
# blas_set_num_threads(1)


source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts_omega_flex.R")
source("code/nuts/adnuts_helper.R")

source("code/hyper_beta_omega.R")
source("code/hyperpar.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

sourceCpp("code/rcpp/posterior_grad_xx2.cpp")

# library(devtools)
# install_github("https://github.com/cran/qrLMM")
library(qrLMM)
data("Cholesterol", package="qrLMM")

data <- Cholesterol

# rescale
data$age <- rescale(Cholesterol$age, to = c(0, 1))
data$year <- rescale(Cholesterol$year, to = c(0, 1))

#standardize
data$cholst <- scale(data$cholst)

seed <- 42
set.seed(seed)
source("code/nuts/nuts_sd_te.R")



# construct centered precision matrices for tensor spline for elicitation of theta
q1 <- q2 <- 10

I <- diag(q2)

K <- vector("list", 2)
P <- diff(diag(q1 - 1), difference = 1)
Pm1 <- matrix(0, q1 - 1, q1)
Pm1[2:(q1 - 1), 2:q1] <- P
K[[1]] <- Pm1 %x% I
K[[1]] <- t(K[[1]])%*%K[[1]]


I2 <- diff(diag(q2), difference = 1)
I21 <- diff(diag(q2), difference = 2)
I1 <- diag(q1)
K[[2]] <- matrix(0, q2 - 2 + (q1 - 1) * (q2 - 1), q1 * q2)
K[[2]][1:(q2 - 2), ] <- t(I1[1, ]) %x% I21
K[[2]][(q2 - 1):nrow(K[[2]]), ] <- I1[2:q1, ] %x% I2
K[[2]] <- t(K[[2]])%*%K[[2]]

Zc <- diag(q1 * q2)
Zc <- Zc[, -q2]
D1 <- t(diff(diag(q2)))
#Zc[1:q2, 1:(q2 - 1)] <- D1

K1 <- t(Zc)%*%K[[1]]%*%Zc
K2 <- t(Zc)%*%K[[2]]%*%Zc
K <- K1 + K2

theta <- hyper_beta_omega(K1, K2, ct=3, alpha=0.01)$root
theta # !!!value of theta is currently hard-coded to this value!!!


object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-6) +
                    hx_lin(sex) + hx_lin(year),
                  family = "gaussian", data=data,
                  iterations = 2000,
                  hyperparams = list(ct=3, alpha=0.1),
                  nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# save(object_te, file="processed_data/fram_te_m10_sd_2000.RData")

# Final acceptance ratio=0.95, and target=0.95
# Final step size=0.013; after 1000 warmup iterations
# Elapsed Time: 479.7 seconds (Warmup)
# Elapsed Time: 431.4 seconds (Sampling)
# Elapsed Time: 911.1 seconds (Total)
# Warning message:
#   In NUTS(n_iter = iterations, warmup = warmup, burnin = burnin, xx = xx,  :
#             max_treedepth(12) reached

load("processed_data/framingham/fram_te_m10_sd_2000.RData")

object_te2$IC
