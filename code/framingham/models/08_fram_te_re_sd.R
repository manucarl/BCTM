##
## Script name: 08_fram_te_re_sd.R
##
## Purpose of script: estimates model bctm_tensor_re_sd for Framingham
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

source("code/hyper_omega_beta.R")

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
# data$age <- scale(data$age)

seed <- 42
set.seed(seed)
source("code/nuts/nuts_sd_te_re.R")



object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-6) +
                     hx_lin(sex) + hx_lin(year) + hx_re(newid, data=data),
                   family = "gaussian", data=data,
                   iterations = 4000,
                   nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# Chain 1, Iteration: 2000/2000 [100%] (Sampling)                                                                                                                          
# Final acceptance ratio=0.91, and target=0.95
# Final step size=0.01; after 1000 warmup iterations
# Elapsed Time: 4611.5 seconds (Warmup)
# Elapsed Time: 6330.8 seconds (Sampling)
# Elapsed Time: 10942.4 seconds (Total)

# save(object_te, file="processed_data/fram_te_sd_m10_re.RData")
load("processed_data/framingham/fram_te_m10_re_sd.RData")

object_te$IC
# $DIC
# DIC       pD     Dbar     Dhat 
# 1568.255  214.066 1354.189 1140.123 
# 
# $WAIC1
# WAIC1     llpd   pwaic1 
# 1516.996 -595.691  162.807 
# 
# $WAIC2
# WAIC2     llpd   pwaic2 
# 1576.949 -595.691  192.784 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate        SE
# elpd_loo -795.1003 27.414879
# p_loo     199.4094  9.315036
# looic    1590.2007 54.829759

