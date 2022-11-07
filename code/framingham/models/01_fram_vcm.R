##
## Script name: 01_fram_vcm.R
##
## Purpose of script: estimates model bctm_vcm for Framingham
##
## Author: BLIND
##
## Date Created: 2022-10-5
##
## Email: BLIND
##
## ---------------------------

library(dplyr)

source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts.R")
source("code/nuts/adnuts_helper.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "qrLMM",
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

##########################################################################################################################################################
# model 1 - varying coefficient for age -----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_vcm <- bctm(cholst ~  hy_sm(cholst, data=data) +
                     hyx_vcm(cholst, by=age, data=data, add_to_diag=10e-6) +
                     hx_lin(sex) + hx_lin(year),
                   family = "gaussian", data=data,
                   iterations = 4000,
                   hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# Final acceptance ratio=0.96, and target=0.95
# Final step size=0.003; after 2000 warmup iterations
# Elapsed Time: 1715.5 seconds (Warmup)
# Elapsed Time: 2218.4 seconds (Sampling)
# Elapsed Time: 3933.9 seconds (Total)


# save(object_vcm, file="processed_data/framingham/fram_vcm.RData")
load("processed_data/framingham/fram_vcm.RData")

# information criteria (we refer to WAIC2 in the paper)
object_vcm$IC
# $DIC
# DIC       pD     Dbar     Dhat 
# 2760.335    9.461 2750.874 2741.413 
# 
# $WAIC1
# WAIC1      llpd    pwaic1 
# 2760.465 -1370.641     9.591 
# 
# $WAIC2
# WAIC2      llpd    pwaic2 
# 2760.811 -1370.641     9.764 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate         SE
# elpd_loo -1380.444796 25.5413439
# p_loo        9.803479  0.7736101
# looic     2760.889593 51.0826879