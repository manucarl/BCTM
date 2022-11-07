##
## Script name: vet_models.R
##
## Purpose of script: calculates contents of Table 1 (VA lung cancer) in supplement
##
## Author: BLIND
##
## Date Created: 2020-10-5
##
## Email: BLIND
##
## ---------------------------

library(dplyr)
rm(list = ls())


source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_fun.R")
source("code/bctm_design.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts.R")
source("code/nuts/adnuts_helper.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "caret",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots", "survival")
load_inst(packages)

sourceCpp("code/rcpp/posterior_grad_xx2.cpp")


data(cancer, package="survival")

data2 <- veteran[!as.logical(veteran$prior),]

dummies <- predict(dummyVars(~ celltype, data = data2), newdata = data2)
colnames(dummies) <- c("squam", "small", "adeno", "large")
data2 <- cbind(data2[,c("time", "status", "karno", "trt")], dummies)

# data2$karno <- rescale(data2$karno, c(0,1))
data2$trt <- ifelse(data2$trt == 2, 1, 0)
# data2$karno <- scale(data2$karno) 

n <- nrow(data)


family <- "logistic"

seed <- 42
# model 1 - po with censoring-------------------------------------------------------------------------------------------------------------------------
object_lin <- bctm(time ~ hy_sm(time,data=data2, add_to_diag=10e-6, center=T) + hx_lin(karno) + hx_lin(adeno) + hx_lin(small) + hx_lin(squam),  
                  family = family, data=data2,
                  cens=as.logical(data2$status),
                  iterations = 4000, 
                  intercept=T,
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)
# save(object_lin, file="processed_data/veteran/vet_po_lin.RData")
load("processed_data/veteran/vet_po_lin.RData")

object_lin$IC
 

burnin <- object_lin$mcmc$burnin

# indices of exponentiated coefficients
exp_ident <- object_lin$model$exp_ident


# beta samples
beta_samples <- object_lin$samples$beta[(2000+1):4000,]

# beta_tilde samples
bt_samples <- beta_samples
bt_samples[,exp_ident] <- exp(beta_samples[,exp_ident])


# same model with MLT
library(mlt)
data3 <- data2
data3$karno <- as.numeric(data3$karno)
data3$time <- with(data2, survival::Surv(time, status))
var_t <- numeric_var("time", support = c(0, 500), bounds = c(0, 600))
b_t <- Bernstein_basis(var_t, order = 10, ui = "increasing")
b_R <- as.basis(~ karno  + adeno + small + squam, data = data3, remove_intercept = TRUE, negative=F)
ctm_t <- ctm(response = b_t, shifting = b_R, todistr = "Logistic")
mlt_t <- mlt(ctm_t, data = data3, scale = TRUE)


coefs <- tibble(Parameter = c("score", "adeno vs. large", "small vs. large", "squamous vs. large"),
                BCTM= round(apply(bt_samples, 2, median), 3)[-(1:20)],
                MLT =round(coef(mlt_t), 3)[-(1:11)],
                MPT = c(-0.055, 1.303, 1.362, -0.173) # values from paper
)

sds <- tibble(BCTM_sd=round(apply(bt_samples, 2, sd), 3)[-(1:20)],
              MLT_sd= round(sqrt(as.vector(diag(vcov(mlt_t)))), 3)[-(1:11)],
              MPT_sd=   c(0.010, 0.559, 0.527, 0.580)
)


table1 <- bind_cols(coefs, sds, )[,c(1,2,5, 3,6, 4,7)]
table1
# A tibble: 4 × 7
# Parameter            BCTM BCTM_sd    MLT MLT_sd    MPT MPT_sd
# <chr>               <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#   1 score              -0.056   0.011 -0.057  0.011 -0.055  0.01 
# 2 adeno vs. large     1.35    0.591  1.36   0.561  1.30   0.559
# 3 small vs. large     1.45    0.551  1.46   0.533  1.36   0.527
# 4 squamous vs. large -0.186   0.626 -0.188  0.598 -0.173  0.58 

xtable::xtable(table1, digits=3)

# model 2 - po with censoring and nonlinear effect for karnofsky score-------------------------------------------------------------------------------------------------------------
object_nl <- bctm(time ~ hy_sm(time,data=data2, q=22, add_to_diag=10e-6, center=T) + hx_sm(karno, data=data2, q=6) + hx_lin(adeno) + hx_lin(small) + hx_lin(squam),  
            family = family, data=data2,
            cens=as.logical(data2$status),
            iterations = 4000,
            intercept=T,
            hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# save(object_nl, file="processed_data/veteran/vet_po_nonlin.RData")
load("processed_data/veteran/vet_po_nonlin.RData")

object_nl$IC

 