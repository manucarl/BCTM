##
## Script name: 02_fram_vcm_sd.R
##
## Purpose of script: estimates model bctm_vcm_sd for Framingham
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

data$cholst <- scale(data$cholst)%>% as.vector


seed <- 123
set.seed(seed)



source("code/nuts/nuts_sd.R")

source("code/hyper_beta.R")
source("code/hyper_omega_beta.R")


object_vcm_sd <- bctm(cholst ~  hy_sm(cholst, data=data, add_to_diag=10e-6) +
                        hyx_vcm(cholst, by=age, center=F, data=data, q=20, add_to_diag=10e-6) +
                        hx_lin(sex) + hx_lin(year),
                      family = "gaussian", data=data,
                      iterations = 4000, warmup=2000, burnin=2000,
                      nuts_settings=list(adapt_delta = 0.99, max_treedepth=12), seed = seed)



# save(object_vcm_sd, file="processed_data/framingham/fram_vcm_sd.RData")
load("processed_data/framingham/fram_vcm_sd.RData")
object_vcm_sd$IC
