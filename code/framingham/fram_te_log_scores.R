##
## Script name: fram_te_log_scores.R
##
## Purpose of script: estimates log scores for Framingham tensor model
##
## Author: BLIND
##
## Date Created: 2022-10-8
##
## Email: BLIND
##
## ---------------------------

library(dplyr)


library(dplyr)

source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts.R")
source("code/nuts/adnuts_helper.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "coda", "qrLMM", "loo",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots", "microbenchmark", "caret")


load_inst(packages)


sourceCpp("code/rcpp/posterior_grad_xx2.cpp")

data("Cholesterol", package="qrLMM")

all_data <- Cholesterol

# scale for better sampling
all_data$age <- rescale(Cholesterol$age, to = c(0, 1))
all_data$year <- rescale(Cholesterol$year, to = c(0, 1))


# function that scales prediction variable according to the scaling in training data
scale_pred <- function(x, var){
  (x - attr(var, "scaled:center"))/
    attr(var, "scaled:scale")
}

seed <- 42
set.seed(seed)

folds <- createFolds(all_data$cholst, k = 10)

##########################################################################################################################################################
# model 2 - tensor spline for age -----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################

log_scores_te <- mclapply(folds, function(holdout) {

data <- all_data[-holdout,]

test_data <- all_data[holdout,]
object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-6) +
                    hx_lin(sex) + hx_lin(year),
                  family = "gaussian", data=data,
                   iterations = 5000, warmup = 1000, burnin = 1000,
                  hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)

# indices of parameters that are exp-transformed
exp_ident <-object_te$model$exp_ident

# beta samples
betas <- object_te$samples$beta

#beta_tilde samples
bts <- betas
bts[,exp_ident] <- exp(bts[,exp_ident])

# prediction df
pred_grid <- test_data

# 
predvars <- object_te$predvars[["hyx_sm(cholst, age)"]]

#  tensor prediction matrix
B_pred <-predvars$B(pred_grid)
Bp_pred <-predvars$Bp(pred_grid)

# design matrices for predictions
X_pred <- cbind(1, B_pred, 0, pred_grid$year)
Xp_pred <- cbind(0, Bp_pred, 0, 0)

# vector of posterior means of reparametrized basis coefficients beta_tilde 
bt <- object_te$beta_tilde

# estimated conditional transformation function
h_est <- X_pred%*%bt


# estimated derivative of conditional transformation function
# hp_est <- Xp_pred%*%bt
# estimated conditional density
# d_est <- dnorm(h_est)*hp_est
 
pnorm(h_est, log.p = TRUE)

}, mc.cores = 10
)

save(log_scores_te, file = "processed_data/framingham/fram_te_log_scores.RData")
