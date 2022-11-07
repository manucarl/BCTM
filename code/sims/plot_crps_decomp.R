##
## Script name: plot_crps_decomp.R
##
## Purpose of script: plots CRPS decomposition for scenario with zero and five noise parameters from BCTM and MLT
##
## Author: BLIND
##
## Date Created: 2022-10-7
##
## Email: BLIND
##
## ---------------------------

library(tidyverse)
library(dplyr)
library(tidyr)


source("code/bctm_fun.R")
source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts_omega.R")
source("code/nuts/adnuts_helper.R")

source("code/mlt_sims/simfuns.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "bamlss",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

pnon <- 0

# e1 <- hyx_sm(y, x1, data=data)
# e2 <- hyx_sm(y, x2, data=data)
# e3 <- hyx_sm(y, x3, data=data)
# e4 <- hyx_sm(y, x4, data=data)
# e5 <- hyx_sm(y, x5, data=data)
# e6 <- hyx_sm(y, x6, data=data)
# e7 <- hyx_sm(y, x7, data=data)
#   
# 
# effects <- list(e1, e2, e3, e4, e5, e6, e7)

Bpreds <- Bppreds <- vector("list", pnon+2)

X_mats <- vector("list", nrow(pred_grid))

# want ~ 1000 samples of estimated cdfs
Bpreds_list <- Bppreds_list <- vector("list", pnon+2)


# for(i in 1:nrow(pred_grid)){
#   xsmat <- as.matrix(pred_grid[i, colnames(pred_grid) %in% varn], ngrid, 1)
#   
#   # data <- expand.grid(ys, xsmat)
#   
#   for(j in 1:(pnon+2)){
#     pred_data <- expand.grid(ys, xsmat[j])
#     colnames(pred_data) <- c("y", varn[j])
#     
#     Bpreds_list[[j]] <- effects[[j]]$predvars$B(pred_data)
#     Bppreds_list[[j]] <- effects[[j]]$predvars$Bp(pred_data)
#     
#     
#   }
#   
#   
#   Xpred <- cbind(1,do.call(cbind, Bpreds_list))
#   
#   X_mats[[i]] <- Xpred
#   
#   
# }

compute_quantile <- function(alpha){
  
  h <- function(q, data){
    
    pnon <- ncol(data)-2
    pred_data <- data.frame(y=q, data)
    for(j in 1:(pnon+2)){
      # pred_data <- test_grid[, c(paste0("x", j), "y")]
      Bpreds_list[[j]] <- effects[[j]]$predvars$B(pred_data)
      
    }
    
    Xpred <- cbind(1,do.call(cbind, Bpreds_list))
    
    # 3) Abweichung von alpha bestimmen
    h_est <- Xpred%*%bt
    h_est
  } 
  
  tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
  
  
  inverse <- function(f, data,lower, upper){
    function(y){
      uniroot(function(x){f(x, data) - y}, lower = lower, upper = upper, tol=1e-3)[1]
    }
  }
  
  q_preds <- rep(0, nrow(pred_grid))
  
  for(i in 1:nrow(pred_grid)){
    h_inverse <- inverse(h, data=pred_grid[i,1:(pnon+2)],-100, 100)
    q_preds[i] <-h_inverse(qnorm(alpha))$root
  }
  
  
  q_preds
}

quantile_levels <- seq(0.05, 0.95, by= 0.01)
#conditional_quantiles_bctm  <- sapply(quantile_levels, function(x) compute_quantile(x)) %>%
#  as_tibble %>% set_names(seq(0.05, 0.95, by= 0.01))

gg <- lapply(c(0, 5), function(p) {
  
  load(paste0("processed_data/sims/score_sims_nonlin_pnon", pnon, "_its4000_burnin2000_m8.RData"))
  
  data <- res[[1]][[1]]$data
  
  varn <- colnames(data)
  varn <- varn[grep("^x", varn)]
  
  
  # construct grid for prediction check
  pred_grid <- dgp(pnon = pnon)
  ngrid <- 100
  ys <- seq(from = min(pred_grid$y), to = max(pred_grid$y), length = ngrid)
  
  
  y<- data$y
  
load(paste0("processed_data/sims/score_sims_nonlin_pnon", p, "_its4000_burnin2000_m8.RData"))
  
  
load(paste0("processed_data/sims/conditional_quantiles_bctm_p", p, "_rep1_bctm.Rdata"))
  
crps_decomp <- function(q, cq_data)  -2*mean(unlist((as.numeric(pred_grid$y < cq_data[,q]) - quantile_levels[q])*(cq_data[,q] - pred_grid$y)))

crps_decomp_bctm  <- sapply(seq_along(quantile_levels), crps_decomp,  conditional_quantiles_bctm)




mlt_object <- res[[1]][[1]]$mlt

conditional_quantiles_mlt <- t(predict(mlt_object, newdata=pred_grid, type="quantile", prob=quantile_levels))

crps_decomp_mlt <- sapply(seq_along(quantile_levels), crps_decomp,  conditional_quantiles_mlt)



# f <- paste0("y ~ s(x1) + s(x2) + s(x3)+ s(x4)+ s(x5)+ s(x6)+ s(x7)"
if (length(varn) > 0) {
  f <- paste("y ~ ", paste(paste0("s(", varn, ")"), 
                           collapse = "+")) %>% as.formula
}

## Quantile regression.
bamlss_models <- lapply(quantile_levels, function(alpha) bamlss(f,
                                                           data = data, optimizer = FALSE, #sampler = BayesX,
                                                           family = gF("quant", prob = alpha))
)


conditional_quantiles_bamlss_qr <- sapply(seq_along(quantile_levels), function(i){
  
  q_alpha <- predict(bamlss_models[[i]], newdata=pred_grid)
  q_alpha 
}) %>% as_tibble %>% set_names(quantile_levels)

crps_decomp_bamlss_qr <- sapply(seq_along(quantile_levels), crps_decomp,  conditional_quantiles_bamlss_qr)

gg <- tibble(Quantile = quantile_levels, bind_rows("bctm" = crps_decomp_bctm, "mlt" = crps_decomp_mlt, "bamlss_qr" = crps_decomp_bamlss_qr, .id="model")) %>% 
  pivot_longer(cols=-Quantile, names_to="Model")

gg
})

# save(gg, file="processed_data/quantile_scores_gg.Rdata")
load("processed_data/sims/quantile_scores_gg.Rdata")

gg <- bind_rows(gg, .id="p")
gg$Model <- factor(gg$Model, levels=c("bctm", "mlt", "bamlss_qr"), labels=c("Full BCTM", "Full MLT", "BAMLSS QR"))
p1 <- gg %>% ggplot() +
  geom_line(aes(x=Quantile, y=-value, colour=Model)) +
  facet_grid(~factor(p, labels=c("p = 0", "p = 5"))) +
  ylab("Quantile score") +
  theme(legend.position="bottom",
        strip.background = element_rect(colour="black"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, size=20),
        axis.text.y = element_text(angle = 45, size=20))
p1

ggsave("manuscript/figs/crps_decomp.png", plot=p1, height=6)
# points(quantile_levels ,sapply(seq_along(quantile_levels), crps_decomp,  conditional_quantiles_bamlss_qr), type="l")

