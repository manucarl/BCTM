##
## Script name: 06_fram_vcm_re_sd.R
##
## Purpose of script: estimates model bctm_vcm_re_sd for Framingham
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
source("code/nuts/nuts_sd.R")
source("code/hyper_beta.R")
source("code/hyper_omega_beta.R")
##########################################################################################################################################################
# model 2 - varying coefficient for age plus random effect -----------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_vcm_re_sd <- bctm(cholst ~  hy_sm(cholst, data=data) +
                        hyx_vcm(cholst, by=age, data=data, add_to_diag=10e-6) +
                        hx_lin(sex) + hx_lin(year) + hx_re(newid, data=data),
                      family = "gaussian", data=data,
                      iterations = 4000,
                      nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)



# save(object_vcm_re_sd, file="processed_data/framingham/fram_vcm_re_sd.RData")

load("processed_data/framingham/fram_vcm_re_sd.RData")

object_vcm_re_sd$IC
# $DIC
# DIC       pD     Dbar     Dhat 
# 1571.726  191.786 1379.940 1188.154 
# 
# $WAIC1
# WAIC1     llpd   pwaic1 
# 1525.581 -617.150  145.641 
# 
# $WAIC2
# WAIC2     llpd   pwaic2 
# 1578.538 -617.150  172.119 
# 
# $WAIC3
# WAIC3     llpd   pwaic3 
# 1234.299 -617.150    0.000 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate        SE
# elpd_loo -795.9944 27.410806
# p_loo     178.8448  8.609026
# looic    1591.9888 54.821613