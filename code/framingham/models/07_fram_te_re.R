##
## Script name: 07_fram_te_re.R
##
## Purpose of script: estimates model bctm_tensor_re for Framingham
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
                    hx_lin(sex) + hx_lin(year) + hx_re(newid, data=data),
                  family = "gaussian", data=data,
                  iterations = 4000,
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# save(object_te, file="processed_data/framingham/fram_te_m10_re.RData")
load("processed_data/framingham/fram_te_m10_re.RData")

object_te$IC
# $DIC
# DIC       pD     Dbar     Dhat 
# 1562.324  192.075 1370.250 1178.175 
# 
# $WAIC1
# WAIC1     llpd   pwaic1 
# 1518.670 -610.915  148.420 
# 
# $WAIC2
# WAIC2     llpd   pwaic2 
# 1573.975 -610.915  176.072 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate       SE
# elpd_loo -793.5353 27.63926
# p_loo     182.6203  8.97875
# looic    1587.0706 55.27852

