##
## Script name: fram_log_scores.R
##
## Purpose of script: performs 10 fold CV and calculates log scores for each fold for both framingham models
##
## Author: Manuel Carlan
##
## Date Created: 2021-09-27
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

source("bctm_utils.R")
source("bctm_design_funs2.R")
source("bctm_design.R")
source("bctm_fun.R")

sourceCpp("rcpp/posterior_grad_xx2.cpp")
library(fastDummies)
library(tidyverse)
# get data
data("LeukSurv", package="spBayesSurv")
data <- LeukSurv[order(LeukSurv$district), ]
n <- nrow(data) 

# data$time <- scale(data$time)
#get boundary file
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))

# construct map in graph format for neighbour matrix
nmat <- bnd2gra(nwengland)


nmat_centroids <- get.centroids(nwengland)

# scale for faster sampling
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)


# data$time <- scale(data$time)
# minimum-extreme-value distribution for proportional hazards
family <- "mev"



seed <- 123

set.seed(42)

# folds <- createFolds(data$time, k = 10)

folds <- fold(data, k=10, cat_col="district") %>% ungroup %>% 
  mutate(fold = as.numeric(as.character(.folds)))



source("nuts/nuts.R")

log_scores_ph <- lapply(1:10, function(i) {
  
  test_data = folds %>% filter(fold == i)
  train_data <<- folds %>% filter(fold != i)
  
  
  
  load(paste0("../processed_data/leukemia/leuk_ph_log_score_batch_", i, ".RData"))
  
  
  # # indices of parameters that are exp-transformed
  exp_ident <-object_ph$model$exp_ident
  
  # beta samples
  #betas <- object_ph$samples$beta
  
  #beta_tilde samples
  #bts <- betas
  #bts[,exp_ident] <- exp(bts[,exp_ident])
  
  # prediction df
  pred_grid <- test_data
  
  X_spat <- dummy_cols(pred_grid, select_columns = "district", remove_selected_columns=T)
  
  X_spat <- X_spat[, (ncol(X_spat)-length(pred_grid$district %>% unique)-1):ncol(X_spat)][,-(1:2)]
  # vector of posterior means of reparametrized basis coefficients beta_tilde
  bt <- object_ph$beta_tilde
  #
  predvars <- object_ph$predvars[["hy_sm(time)"]]
  
  #  tensor prediction matrix
  B_pred <-predvars$B(pred_grid)
  Bp_pred <-predvars$Bp(pred_grid)
  
  # design matrices for predictions
  X_pred <- cbind(B_pred, X_spat, pred_grid$sex, pred_grid$wbc, pred_grid$tpi) %>% as.matrix
  #Xp_pred <- cbind(0, Bp_pred, 0, 0)
  

  
  # estimated conditional transformation function
  h_est <- X_pred%*%bt
  
  
  # estimated derivative of conditional transformation function
  # hp_est <- Xp_pred%*%bt
  # estimated conditional density
  # d_est <- dnorm(h_est)*hp_est
  
  pnorm(h_est, log.p = TRUE)
  
})

object_ph$
log_scores_ph %>% unlist %>% sum


library(spBayesSurv)

nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))


mcmc <- list(nburn = 2000, nsave = 2000, nskip = 1, ndisplay = 1000)
prior <- list(maxL = 40)
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)
# Enew <- E %>% as.matrix
# class(Enew) <- "matrix"
res4 <- survregbayes(formula = Surv(time) ~ +age + sex+ wbc +tpi +
                       frailtyprior("car", district), data = data, survmodel = "PH",
                     dist = "loglogistic", mcmc = mcmc, Proximity = E, scale.designX=FALSE, InitParamMCMC=FALSE)

res4$WAIC
log_scores_te_sd <- lapply(seq_along(folds), function(i) {
  
  test_data = folds %>% filter(fold == i)
  train_data <<- folds %>% filter(fold != i)
  
  
  
  load(paste0("../processed_data/leukemia/leuk_ph_bayessurv_log_score_batch_", i, ".RData"))
  
  plot(res, xnewdata=test_data[1,])
  
  # estimated derivative of conditional transformation function
  # hp_est <- Xp_pred%*%bt
  # estimated conditional density
  # d_est <- dnorm(h_est)*hp_est
  
  pnorm(h_est, log.p = TRUE)
  
})




library(dplyr)


load("processed_data/framingham/fram_vcm.RData")
#object_vcm$IC$WAIC1

load("processed_data/framingham/fram_vcm_log_scores.RData")
#log_scores_vcm %>% unlist %>%  sum


VCM <- c(Model = "VCM", object_vcm$IC$WAIC1, log_scores = log_scores_vcm %>% unlist %>%  sum)


load("processed_data/framingham/fram_vcm_re.RData")
#object_vcm_re$IC$WAIC1



load("processed_data/framingham/fram_vcm_log_scores_sd_c3.RData")
#log_scores_vcm %>% unlist %>%  sum





load("processed_data/framingham/fram_te.RData")
#object_te$IC$WAIC1

#log_scores_te %>% unlist %>% sum

Tensor <- c(Model = "Tensor", object_te$IC$WAIC1, log_scores = log_scores_te %>% unlist %>% sum)


load("processed_data/framingham/fram_te_sd.RData")
#object_te$IC$WAIC1

#log_scores_te_sd %>% unlist %>% sum

Tensor_sd <- c(Model = "Tensor_sd", object_te$IC$WAIC1, log_scores = log_scores_te_sd %>% unlist %>% sum)


rbind(VCM, Tensor, Tensor_sd)
load("processed_data/framingham/fram_te_re.RData")
object_te_re$IC$WAIC1

load("processed_data/framingham/fram_te_re_sd.RData")
object_te$IC$DIC








# 
# # sd prior
# load("processed_data/fram_vcm_log_scores_sd_c1.RData")
# log_scores_vcm %>% unlist %>%  sum
# 
# 
# 
# load("processed_data/fram_vcm_log_scores_sd_c3_alpha01.RData")
# log_scores_vcm %>% unlist %>%  sum
# 
# 
# load("processed_data/fram_te_log_scores_sd_c3_alpha005.RData")
# log_scores_te %>% unlist %>%  sum
# 
# 
# load("processed_data/framingham/fram_te_log_scores.RData")
# log_scores_te %>% unlist %>%  sum
# 
# load("processed_data/framingham/fram_vcm_log_scores.RData")
# log_scores_vcm%>% unlist %>%  sum
