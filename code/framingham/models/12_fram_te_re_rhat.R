
library(dplyr)
# rm(list = ls())
path <- "/home/mcarlan/bctm_jasa_revision/code"
setwd(path)

source("bctm_utils.R")
source("bctm_design_funs2.R")
source("bctm_design.R")
source("bctm_fun.R")

source("nuts/nuts_utils.R")
source("nuts/adnuts_helper.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "coda",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots", "microbenchmark", "caret")
load_inst(packages)


library("RhpcBLASctl")
omp_set_num_threads(1)
blas_set_num_threads(1)


sourceCpp("rcpp/posterior_grad_xx2.cpp")

data("Cholesterol", package="qrLMM")

data <- Cholesterol

# scale for better sampling
data$age <- rescale(Cholesterol$age, to = c(0, 1))
data$year <- rescale(Cholesterol$year, to = c(0, 1))
data$cholst <- scale(data$cholst)

seed <- 123

# function that scales prediction variable according to the scaling in training data
scale_pred <- function(x, var){
  (x - attr(var, "scaled:center"))/
    attr(var, "scaled:scale")
}

source("nuts/nuts.R")
#source("hyper_omega_beta.R")

##########################################################################################################################################################
# model 2 - tensor spline for age -----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_te <- mclapply(1:4, function(x) bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-6) +
                    hx_lin(sex) + hx_lin(year)+ hx_re(newid, data=data),
                  family = "gaussian", data=data, start = "random",
                  iterations = 4000, warmup = 2000, burnin = 2000, 
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.99, max_treedepth=12), seed = seed),
mc.cores = 4
)

save(object_te, file = "../fram_te_rhat_runs.RData")