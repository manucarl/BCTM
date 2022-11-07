##
## Script name: 05_fram_vcm_re.R
##
## Purpose of script: estimates model bctm_vcm_re for Framingham
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

##########################################################################################################################################################
# model 2 - varying coefficient for age plus random effect -----------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_vcm_re <- bctm(cholst ~  hy_sm(cholst, data=data) +
                        hyx_vcm(cholst, by=age, data=data, add_to_diag=10e-6) +
                        hx_lin(sex) + hx_lin(year) + hx_re(newid, data=data),
                      family = "gaussian", data=data,
                      iterations = 4000,
                      hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# Final acceptance ratio=0.91, and target=0.95
# Final step size=0.005; after 2000 warmup iterations
# Elapsed Time: 19893.1 seconds (Warmup)
# Elapsed Time: 28720.6 seconds (Sampling)
# Elapsed Time: 48613.6 seconds (Total)


# save(object_vcm_re, file="processed_data/framingham/fram_vcm_re.RData")
load("processed_data/framingham/fram_vcm_re.RData")

object_vcm_re$IC
# $DIC
# DIC       pD     Dbar     Dhat 
# 1569.356  191.505 1377.852 1186.347 
# 
# $WAIC1
# WAIC1     llpd   pwaic1 
# 1524.039 -615.832  146.188 
# 
# $WAIC2
# WAIC2     llpd   pwaic2 
# 1577.198 -615.832  172.767 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate        SE
# elpd_loo -795.2175 27.676670
# p_loo     179.3856  8.586777
# looic    1590.4350 55.353339