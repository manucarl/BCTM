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
path <- "D:/onedrive/bctm_jasa_revision/code/"
setwd(path)

source("bctm_utils.R")
source("bctm_design_funs2.R")
source("bctm_design.R")
source("bctm_fun.R")

source("nuts/nuts_utils.R")
source("nuts/nuts.R")
source("nuts/adnuts_helper.R")

packages <- c("spBayesSurv", "survival", "BayesX", "scam", "Matrix", "Rcpp", "RcppArmadillo", "MCMCpack", "sf", "rgeos", "ggplot2")
load_inst(packages)


sourceCpp("rcpp/posterior_grad_xx2.cpp")

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
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)
data$tpi <- scale(data$tpi)

# data$sex[data$sex ==0 ] <- -
data$sex[data$sex ==0 ] <- -1

# data$time <- scale(data$time)
# minimum-extreme-value distribution for proportional hazards
family <- "gaussian"


seed <- 1
its <- 2000



mcmcplots::mcmcplot(object_ph$samples$beta)


load("../processed_data/leukemia/leuk_ph.RData")
object_ph$IC

source("nuts/nuts.R")

object_ph0 <- bctm(time ~ hy_sm(time,data=data,  center=T, q = 40, add_to_diag=10e-6)+
                     # hx_spat(district, data=data, nmat=nmat)+
                    hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) , #cens = as.logical(data$cens),
                  family = "mev", data=data, iterations = 2000, intercept=T, # remove intercept 
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

object_ph0$IC
# proportional hazards model with spatial effect -----------------------------------------------------------------------------------------------
object_ph <- bctm(time ~ hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                    hx_spat2(district, data=data, nmat=nmat)+
                    hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi), #, cens = as.logical(data$cens),
                  #start=rep(1, 47),
                  family = family, data=data, iterations = 2000, intercept=T, # remove intercept , star
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)
object_ph$IC
bt_samples <-  object_ph$samples$beta

exp_ident <- object_ph$model$exp_ident
bt_samples[,exp_ident] <- exp(bt_samples[,exp_ident])

beta_tilde <- colMeans(bt_samples)

X <- object_ph$X
Xp <- object_ph$Xp

h <- X%*%t(bt_samples)
hp <- Xp%*%t(bt_samples)

dens <- bctm_dmev
LaplacesDemon::WAIC(dens(h, log=T) + log(hp))

object_ph$IC
mcmcplots::mcmcplot(object_ph$samples$beta)

library(splines)
source("nuts/nuts_omega_spat.R")
# non-proportional hazards model with spatial effect -----------------------------------------------------------------------------------------------
object_nph <- bctm(time ~ hyx_sm(time, age, data=data, q=c(10, 10), center=T)+
                    hx_spat2(district, data=data, nmat=nmat)+
                    hx_lin(sex) + hx_lin(wbc) + hx_lin(tpi) , #cens = as.logical(data$cens),
                  family = family, data=data, iterations = 100, intercept=T, # remove intercept 
                  hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.80, max_treedepth=10), seed = seed)


object_nph$IC
# save(object_nph, file="processed_data/leuk_nph.RData")
load("../processed_data/leukemia/leuk_nph.RData")


object_nph$IC


load("../processed_data/leuk_bayessurv_ph_weibull.RData")



nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))


mcmc <- list(nburn = 2000, nsave = 2000, nskip = 1, ndisplay = 1000)
prior <- list(maxL = 40)
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)
# Enew <- E %>% as.matrix
# class(Enew) <- "matrix"



res4 <- survregbayes(formula = Surv(time) ~ +age + sex+ wbc +tpi,data = data, survmodel = "PH",
                     dist = "lognormal", mcmc = mcmc, Proximity = E, scale.designX=FALSE, InitParamMCMC=FALSE)


res2 <- survregbayes(formula = Surv(time) ~ +age + sex+ wbc +tpi, data = data, survmodel = "PH",
                     dist = "loglogistic", mcmc = mcmc, Proximity = E, scale.designX=FALSE, InitParamMCMC=FALSE)
res3 <- survregbayes(formula = Surv(time) ~ +age + sex+ wbc +tpi, data = data, survmodel = "PO",
                     dist = "loglogistic", mcmc = mcmc, Proximity = E, scale.designX=FALSE, InitParamMCMC=FALSE)


res4 <- survregbayes(formula = Surv(time, cens) ~ +age + sex+ wbc +tpi +
                       frailtyprior("car", district), data = data, survmodel = "PO",
                     dist = "loglogistic", mcmc = mcmc, Proximity = E, scale.designX=FALSE, InitParamMCMC=FALSE)
res4$DIC

res4$WAIC

res5 <- survregbayes(formula = Surv(time) ~ +age + sex+ wbc +tpi +
                       frailtyprior("car", district), data = data, survmodel = "PH",
                     dist = "loglogistic", mcmc = mcmc, Proximity = E, scale.designX=T, InitParamMCMC=FALSE)


save(res2, file="../processed_data/leuk_bayessurv_ph_weibull.RData")


library(mlt)

var_y <- numeric_var("time", support=range(data$time))
B_y <- Bernstein_basis(var_y, order = 20, ui = "increasing")

B_x <- as.basis(~ age + sex + wbc + tpi, data=LeukSurv, remove_intercept=TRUE)

ctm_yx <- ctm(B_y, shifting=B_x,
    data = LeukSurv, todistr = "MinExtrVal", sumconstr=F)

mlt_yx <- mlt(ctm_yx, data = LeukSurv)

B_xs <- var_x <- vector("list", pnon+2)
xvars <- colnames(data)
xvars <- xvars[grep("^x", xvars)]

for(j in 1: (pnon+2)){
  var_x[[j]]<- numeric_var(xvars[j], support = range(data[,xvars[j]]), bounds = range(data[,xvars[j]]))
  B_xs[[j]] <- Bernstein_basis(var_x[[j]], order = m+2, ui = "none")
}
names(B_xs) <- paste0("b", 1: length(B_xs))

ctm_yx <- switch(pnon +1,
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]]),
                     data = data, todistr = "Normal", sumconstr=F),
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]], b3 = B_xs[[3]]),
                     data = data, todistr = "Normal", sumconstr=F),
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]], b3 = B_xs[[3]], b4 = B_xs[[4]]),
                     data = data, todistr = "Normal", sumconstr=F),
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]], b3 = B_xs[[3]], b4 = B_xs[[4]], b5 = B_xs[[5]]), 
                     data = data, todistr = "Normal", sumconstr=F),
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]], b3 = B_xs[[3]], b4 = B_xs[[4]], b5 = B_xs[[5]], b6 = B_xs[[6]]), 
                     data = data, todistr = "Normal", sumconstr=F),
                 ctm(B_y, interacting = c(b1 = B_xs[[1]], b2 = B_xs[[2]], b3 = B_xs[[3]], b4 = B_xs[[4]], b5 = B_xs[[5]], b6 = B_xs[[6]], b7 = B_xs[[7]]), 
                     data = data, todistr = "Normal", sumconstr=F))      


mlt_yx <- mlt(ctm_yx, data = data)

COVMODEL <- ExponentialCovFct()

data$time2 <- gencens(data$time, data$time)
mod <- survspat(formula = Surv(time) ~ age + sex+ wbc +tpi, data = data, dist="weibullHaz",
                shape= get.centroids(nwengland),
                # dist = DIST,
                 mcmc.control = mcmcpars(nits = 500000, burn = 10000, thin = 490)#,
)
                 priors = priors)


y <- rnorm(n, mean=ff1(xm[,1]) + ff2(xm[,2]) + ff3(xm[,3]) + ff4(xm[,4]))
data <- data.frame(y, xm )
colnames(data) <- c("y", paste0("x", 1:4))

library(mlt)
var_y <- numeric_var("y", support = range(data$y), add = c(0, 0),
                     bounds = c(0, Inf))
var_x1 <- numeric_var("x1", support = c(-2, 2), add = c(0, 0),
                      bounds = c(-Inf, Inf))
var_x2 <- numeric_var("x2", support = c(-2, 2), add = c(0, 0),
                      bounds = c(-Inf, Inf))
var_x3 <- numeric_var("x3", support = c(-2, 2), add = c(0, 0),
                      bounds = c(-Inf, Inf))
var_x4 <- numeric_var("x4", support = c(-2, 2), add = c(0, 0),
                      bounds = c(-Inf, Inf))


# stay in Gaussian location-scale trafo by setting order = 1 (linearity)
B_y <- Bernstein_basis(var_y, order =8, ui = "increasing")


B_x1 <- Bernstein_basis(var_x1, order =8, ui = "none")
B_x2 <- Bernstein_basis(var_x2, order = 8, ui = "none")
B_x3 <- Bernstein_basis(var_x3, order = 8, ui = "none")
B_x4 <- Bernstein_basis(var_x4, order = 8, ui = "none")

# set up and estimate mlt
ctm_yx <- ctm(response = B_y, interacting = c(b1= B_x1), shifting= c(b3 = B_x3), todistr = "Normal")
mlt_yx <- mlt(ctm_yx, data = data, scale = T)
