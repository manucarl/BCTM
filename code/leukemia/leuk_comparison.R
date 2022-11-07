## ---------------------------
##
## Script name: leuk01_plots.R
##
## Purpose of script: creates survivor function and spatial plot for leukemia data in "Bayesian Conditional Transformation Models"
##
## Author: BLIND
##
## Date Created: 2020-6-5
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

packages <- c("spBayesSurv", "survival", "BayesX", "scam", "Matrix", "Rcpp", "RcppArmadillo", "MCMCpack", "sf", "rgeos", "ggplot2","RhpcBLASctl", "loo")
load_inst(packages)
# 
# omp_set_num_threads(1)
# blas_set_num_threads(1)

sourceCpp("code/rcpp/posterior_grad_xx.cpp")

# get data
data("LeukSurv", package="spBayesSurv")
data <- LeukSurv[order(LeukSurv$district), ]
n <- nrow(data) 

#get boundary file
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))

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

seed <- 1
its <- 4000


#######################################################################################################
# PH model - bctm
#######################################################################################################

object_ph_cens <- bctm(time ~  hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                          hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi), #+
                        cens = as.logical(data$cens),
                        family = "mev", data=data, iterations = its, intercept=T,
                        hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)


# save(object_ph_cens, file="processed_data/leukemia/leuk_ph_cens.RData")
load("processed_data/leukemia/leuk_ph_cens.RData")
object_ph_cens$IC$WAIC1

#######################################################################################################
# PH model - mlt
#######################################################################################################
library(tram)
AIC(Coxph(Surv(time, event=cens) ~ age + sex+ wbc +tpi , data = data, order=5))
AIC(Coxph(Surv(time, event=cens) ~ age + sex+ wbc +tpi , data = data, order=10))
AIC(Coxph(Surv(time, event=cens) ~ age + sex+ wbc +tpi , data = data, order=20))
AIC(Coxph(Surv(time, event=cens) ~ age + sex+ wbc +tpi , data = data, order=30))


#######################################################################################################
# PH model with random effect - bctm
#######################################################################################################

object_ph_re <- bctm(time ~hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                       hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                       hx_re(district, data),
                     cens = as.logical(data$cens),
                     family = "mev", data=data, iterations = its, intercept=T, 
                     hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)


# save(object_ph_re, file="processed_data/leukemia/leuk_ph_re_cens.RData")
load("processed_data/leukemia/leuk_ph_re_cens.RData")

object_ph_re$IC

#######################################################################################################
# PH model with spatial effect - bctm
#######################################################################################################


object_ph_spat <- bctm(time ~hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                         hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                         hx_spat(district, nmat=nmat, data),
                       cens = as.logical(data$cens),
                       family = "mev", data=data, iterations = its, intercept=F, 
                       hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# save(object_ph_spat, file="processed_data/leukemia/leuk_ph_spat_cens.RData")
load("processed_data/leukemia/leuk_ph_spat_cens.RData")
object_ph_spat$IC$WAIC1$estimates



#######################################################################################################
# NPH model -bctm
#######################################################################################################
# library(splineDesign)

source("code/nuts/nuts_omega.R")

object_nph <- bctm(time ~ hyx_sm(time, age, data=data,  center=T, q = c(10,10), add_to_diag=10e-6)+
                     hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) ,
                   cens = as.logical(data$cens),
                   family = "mev", data=data, iterations = its, intercept=T, 
                   hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)

# save(object_nph, file="processed_data/leukemia/leuk_nph_cens.RData")
load("processed_data/leukemia/leuk_nph_cens.RData")
object_nph$IC$WAIC1$estimates

#######################################################################################################
# NPH model with spatial effect - bctm
#######################################################################################################

source("code/nuts/nuts_omega_spat.R")
object_nph_spat <- bctm(time ~ hyx_sm(time, age, data=data,  center=T, q = c(10,10), add_to_diag=10e-6)+
                          hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                          hx_spat(district, nmat=nmat, data),
                        cens = as.logical(data$cens),
                        family = "mev", data=data, iterations = its, intercept=F, 
                        hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)


# save(object_nph_spat, file="processed_data/leukemia/leuk_nph_spat_cens.RData")
load("processed_data/leukemia/leuk_nph_spat_cens.RData")
object_nph_spat$IC$WAIC1$estimates

