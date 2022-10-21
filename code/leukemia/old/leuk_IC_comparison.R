## ---------------------------
##
## Script name: leuk01_plots.R
##
## Purpose of script: creates survivor function and spatial plot for leukemia data in "Bayesian Conditional Transformation Models"
##
## Author: Manuel Carlan
##
## Date Created: 2020-6-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------
library(dplyr)
# rm(list = ls())
path <- "D:/onedrive/bctm_jasa_revision/"
setwd(path)

source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts.R")
source("code/nuts/adnuts_helper.R")

packages <- c("spBayesSurv", "survival", "BayesX", "scam", "Matrix", "Rcpp", "RcppArmadillo", "MCMCpack", "sf", "rgeos", "ggplot2")
load_inst(packages)


sourceCpp("code/rcpp/posterior_grad_xx.cpp")

# get data
data("LeukSurv", package="spBayesSurv")
data <- LeukSurv[order(LeukSurv$district), ]
n <- nrow(data) 

# data$time <- scale(data$time)
#get boundary file
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))
data("LeukSurv", package="spBayesSurv")

# construct map in graph format for neighbour matrix
nmat <- bnd2gra(nwengland)

# data <- data[data$cens == 0,]

nmat_centroids <- get.centroids(nwengland)

# scale for faster sampling
# data$age <- scale(data$age)%>% as.vector
data$wbc <- scale(data$wbc)%>% as.vector
data$tpi <- scale(data$tpi) %>% as.vector
data$age <- scale(data$age) %>% as.vector
data$sex[data$sex==0] <- -1

# data$time <- scale(data$time)
# minimum-extreme-value distribution for proportional hazards
get_info <- function(object){
  
  X <- object$X
  Xp <- object$Xp
  cind <- data$cens %>% as.logical
  X2 <- X[!cind,]
  X1 <- X[cind,]
  
  Xp2 <- Xp[!cind,]
  Xp1 <- Xp[cind,]
  k <- ncol(X)
  
  xx <- list(exp_ident=object$model$exp_ident-1, X1=X1, X2=X2, Xp1=Xp1)
  
  
  AIC_bctm <- -2*ll_mev_cens(object$posterior_means$beta, xx) + 2*k
  BIC_bctm <- -2*ll_mev_cens(object$posterior_means$beta, xx) + 2*k*log(nrow(X))
  AIC_bctm
  
  BIC_bctm
  
  cat(paste("WAIC: ", object$IC$WAIC1$estimates["waic",1]), "\n",
      paste("AIC: ", AIC_bctm)
  )
  

}



seed <- 1
its <- 2000





object_ph <- bctm(time ~  hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                          hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi), #+
                        cens = as.logical(data$cens),
                        family = "mev", data=data, iterations = 2000, intercept=T,
                        hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

load("processed_data/leukemia/leuk_ph.RData")


get_info(object_ph_cens)
object_ph_cens$beta_tilde

#######################################################################################################
# PH model - mlt
#######################################################################################################
library(tram)
dcmp <- Coxph(Surv(time,cens) ~ age + sex+ wbc +tpi , data = data, order=20)
AIC_mlt <- AIC(dcmp)
AIC_mlt
BIC_mlt <- BIC(dcmp)
BIC_mlt

coef(dcmp)


# fm_data <- Surv(timec, cens) ~ age+  sex+ wbc +tpi 
# B_datay <- Bernstein_basis(var = datay, order = 20, ui = "increasing")
# 
# B_shift <- as.basis(fm_data[-2L], data, remove_intercept=F)
# ph_ctm <- ctm(B_datay, shifting = B_shift, data = data, todistr = "MinExtrVal")
# 
# ph_mlt <- mlt(ph_ctm, data = data, scale=F)
# coef(ph_mlt)
# AIC(ph_mlt)






#######################################################################################################
# PH model with random effect - bctm
#######################################################################################################

object_ph_re <- bctm(time ~hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                       hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                       hx_re(district, data),
                     cens = as.logical(data$cens),
                     family = "mev", data=data, iterations = 2000, intercept=T, 
                     hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

load("processed_data/leukemia/leuk_ph_re.RData")
get_info(object_ph_re)

#object_ph_re$beta_tilde


#######################################################################################################
# PH model with RE -  mlt
#######################################################################################################

library("tramME")
dcmp <- CoxphME(Surv(time,cens) ~ age + sex+ wbc +tpi + (1|district), data = data, order=20)
AIC_mlt <- AIC(dcmp)
AIC_mlt
# BIC_mlt <- BIC(dcmp)
# BIC_mlt

coef(dcmp)


#######################################################################################################
# PH model with spatial effect - bctm
#######################################################################################################

object_ph_spat <- bctm(time ~hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                         hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                         hx_spat2(district, nmat=nmat, data),
                       cens = as.logical(data$cens),
                       family = "mev", data=data, iterations = 2000, intercept=T, 
                       hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

load("processed_data/leukemia/leuk_ph_spat.RData")
get_info(object_ph_spat)

samples <- object_ph_spat$samples$beta %>% as_tibble %>% dplyr::select(starts_with("hx_lin")) %>% slice(1001:2000)
colMeans(samples)
apply(samples, 2, sd)
apply(samples, 2, quantile, c(0.025, 0.975))
apply(samples, 2,summary)

# object_ph_spat$beta_tilde

#######################################################################################################
# NPH model for age - bctm
#######################################################################################################
library(splineDesign)
object_nph <- bctm(time ~ hyx_sm(time, age, data=data,  center=T, q = c(10,10), add_to_diag=10e-6)+
                     hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) ,
                   cens = as.logical(data$cens),
                   family = "mev", data=data, iterations = 2000, intercept=T, 
                   hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=8), seed = seed)


load("processed_data/leukemia/leuk_nph.RData")
get_info(object_nph)

object_nph$beta_tilde

#######################################################################################################
# NPH model for age - mlt
#######################################################################################################

datay <- numeric_var("timec", support = c(1, max(data$time)), bounds = c(0, Inf))
age_var <- numeric_var("age")

data$timec <- with(data, Surv(time , cens))


B_datay <- Bernstein_basis(var = datay, order = 9, ui = "increasing")
B_age <- Bernstein_basis(var = age_var, order=9, ui="none")
fm_data <- Surv(timec, cens) ~   sex+ wbc +tpi 

B_shift <- as.basis(fm_data[-2L], data, remove_intercept=F)
nph_ctm <- ctm(B_datay, shifting = B_shift, interacting = B_age, data = data, todistr = "MinExtrVal")

nph_mlt <- mlt(nph_ctm, data = data, scale=F)
coef(nph_mlt)
AIC(nph_mlt)


#######################################################################################################
# NPH model for age wih spatial effect - bctm
#######################################################################################################

object_nph_re <- bctm(time ~ hyx_sm(time, age, data=data,  center=T, q = c(10,10), add_to_diag=10e-6)+
                          hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                          hx_re(district, data),
                        cens = as.logical(data$cens),
                        family = "mev", data=data, iterations = 2000, intercept=T, 
                        hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

load("processed_data/leukemia/leuk_nph_re.RData")
get_info(object_nph_re)
object_nph_re$beta_tilde

#######################################################################################################
# NPH model for age wih spatial effect - bctm
#######################################################################################################

object_nph_spat <- bctm(time ~ hyx_sm(time, age, data=data,  center=T, q = c(10,10), add_to_diag=10e-6)+
                          hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                          hx_spat2(district, nmat=nmat, data),
                        cens = as.logical(data$cens),
                        family = "mev", data=data, iterations = 2000, intercept=T, 
                        hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

load("processed_data/leukemia/leuk_nph_spat.RData")
get_info(object_nph_spat)
object_nph_spat$beta_tilde




