##
## Script name: fram_log_scores.R
##
## Purpose of script: performs 10 fold CV and calculates log scores for each fold for both framingham models
##
## Author: BLIND
##
## Date Created: 2021-09-27
##
## Email: BLIND
##
## ---------------------------



library(tidyverse)




load("processed_data/framingham/log_scores/fram_vcm_log_scores.RData")
log_scores_vcm %>% unlist %>%  sum

load("processed_data/framingham/log_scores/fram_te_log_scores.RData")
log_scores_te %>% unlist %>%  sum





# sd prior
# load("processed_data/framingham/log_scores/fram_vcm_log_scores_sd_c1.RData")
# log_scores_vcm %>% unlist %>%  sum

load("processed_data/framingham/log_scores/fram_vcm_log_scores_sd_c3.RData")
log_scores_vcm %>% unlist %>%  sum

# load("processed_data/framingham/log_scores/fram_vcm_log_scores_sd_c3_alpha01.RData")
# log_scores_vcm %>% unlist %>%  sum


load("processed_data/framingham/log_scores/fram_te_log_scores_sd_c3_alpha005.RData")
log_scores_te %>% unlist %>%  sum

