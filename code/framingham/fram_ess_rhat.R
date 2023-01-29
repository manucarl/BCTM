##
## Script name: 01_fram_vcm.R
##
## Purpose of script: estimates Rhat and effective sample sizes for Framingham model
##
## Author: Manuel Carlan
##
## Date Created: 2022-10-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

library(dplyr)
library(coda)
library(posterior)


##########################################################################################################
# get Rhat
##########################################################################################################
load("processed_data/framingham/rhat_runs/fram_te_rhat_runs.RData")
mcmc_objects <- lapply(lapply(lapply(object_te, "[[", "samples"), "[[", "beta"), as.mcmc)
combinedchains  <- mcmc.list(mcmc_objects[[1]], mcmc_objects[[2]], mcmc_objects[[3]], mcmc_objects[[4]])

d <- as_draws_rvars(combinedchains)
rhats <- sapply(d, rhat)
all(rhats < 1.1)
mean(rhats)
var(rhats)


load("processed_data/framingham/rhat_runs/fram_te_re_rhat_runs.RData")
mcmc_objects <- lapply(lapply(lapply(object_te, "[[", "samples"), "[[", "beta"), as.mcmc)
combinedchains  <- mcmc.list(mcmc_objects[[1]], mcmc_objects[[2]], mcmc_objects[[3]], mcmc_objects[[4]])

d <- as_draws_rvars(combinedchains)
rhats <- sapply(d, rhat)
all(rhats < 1.1)
mean(rhats)
var(rhats)


load("processed_data/framingham/rhat_runs/fram_te_sd_rhat_runs.RData")
mcmc_objects <- lapply(lapply(lapply(object_te, "[[", "samples"), "[[", "beta"), as.mcmc)
combinedchains  <- mcmc.list(mcmc_objects[[1]], mcmc_objects[[2]], mcmc_objects[[3]], mcmc_objects[[4]])

d <- as_draws_rvars(combinedchains)
rhats <- sapply(d, rhat)
all(rhats < 1.1)
mean(rhats)
var(rhats)


load("processed_data/framingham/rhat_runs/fram_te_re_sd_rhat_runs.RData")
mcmc_objects <- lapply(lapply(lapply(object_te, "[[", "samples"), "[[", "beta"), as.mcmc)
combinedchains  <- mcmc.list(mcmc_objects[[1]], mcmc_objects[[2]], mcmc_objects[[3]], mcmc_objects[[4]])

d <- as_draws_rvars(combinedchains)
rhats <- sapply(d, rhat)
all(rhats < 1.1)
mean(rhats)
var(rhats)


##########################################################################################################
# get ESS
##########################################################################################################


get_ess <- function(object) effectiveSize(object$samples$beta[2001:4000, 1:99] %>% as.mcmc)

# regular tensor model
load("processed_data/framingham/fram_te_m10.RData")

print(paste0("Elapsed time for 1000 samples: ", object_te$samples$time.total / 4 / 60, "minutes"))
get_ess(object_te) %>% mean

# tensor model with scale-dependent prior on tau2
load("processed_data/framingham/fram_te_m10_sd.RData")
print(paste0("Elapsed time for 1000 samples: ", object_te$samples$time.total / 4 / 60, " minutes"))

get_ess(object_te) %>% mean

# regular tensor model with random effect
load("processed_data/framingham/fram_te_m10_re.RData")
print(paste0("Elapsed time for 1000 samples: ", object_te$samples$time.total / 4 / 60, " minutes"))
get_ess(object_te) %>% mean

# tensor model with scale dependent prior an random effect
load("processed_data/framingham/fram_te_m10_re_sd.RData")
print(paste0("Elapsed time for 1000 samples: ", object_te$samples$time.total / 4 / 60, " minutes"))
get_ess(object_te) %>% mean

