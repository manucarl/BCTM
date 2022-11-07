##
## Script name: 03_fram_te.R
##
## Purpose of script: estimates model bctm_tensor for Framingham 
##
## Author: BLIND
##
## Date Created: 2020-10-7
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
source("code/nuts/nuts_omega.R")
source("code/nuts/adnuts_helper.R")

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


object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-6) +
                    hx_lin(sex) + hx_lin(year),
                  family = "gaussian", data=data,
                  iterations = 2000,
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)
# Final acceptance ratio=0.96, and target=0.95
# Final step size=0.012; after 1000 warmup iterations
# Elapsed Time: 424.5 seconds (Warmup)
# Elapsed Time: 544.1 seconds (Sampling)
# Elapsed Time: 968.6 seconds (Total)

# save(object_te, file="processed_data/fram_te_m10.RData")
load("processed_data/framingham/fram_te_m10.RData")

object_te$IC


                        