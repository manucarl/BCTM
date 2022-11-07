##
## Script name: generate_sims_lin.R
##
## Purpose of script: implements linear simulation scenario as shown in the paper
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------

library(dplyr)
seedrepl <- 123

rm(list = ls())
packages <- c("stringr","MASS","doParallel", "loo",  "BayesX","microbenchmark", "Rcpp", "RcppArmadillo", "splines", "mgcv", "Matrix", "MCMCpack", "sdPrior", "R2BayesX",
              "RhpcBLASctl", "scam", "bamlss", "mlt", "rlang", "scoringutils", "scales")


source("code/helpers.R")

load_inst_packages(packages)

# restrict parallel computing 
# library("RhpcBLASctl")
# omp_set_num_threads(1)
# blas_set_num_threads(1)

# code provided by hothorn et al. (2017)
source("code/mlt_sims/simfuns.R")

# posteriors and gradients in rcpp
sourceCpp("code/rcpp/posterior_grad_xx2.cpp")



source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts_lin.R")
source("code/nuts/adnuts_helper.R")



run_sims <- function(seed, its, burnin,  pnon, m, a = 1, b = 0.001){
  
  set.seed(seed)
  
  data <-  dgp(c(100, 100), pnon = pnon)
  n <- nrow(data)
  
  y <- data$y
  
  varn <- colnames(data)
  varn <- varn[grep("^x", varn)]
  
  # select matrix of covariates
  xmat <- xlin <- data %>% dplyr::select(all_of(varn)) %>% as.matrix()
  ymat <- as.matrix( y)
  
  f <- paste("y", "~", "hy_lin(y) +", paste0("hx_lin(" , varn,")", 
                                           collapse=" + "), "+", paste("hyx_lin(y," , varn,")",
                             collapse=" + "))
  mod <- paste("~", "y +", paste0(varn,collapse=" + "), "+", paste("y:" , varn, collapse=" + ")) %>% as.formula
  
  data_train <- data.frame(y, apply(data[,varn], 2, rescale))
  
  res <- bctm(      
    as.formula(f),
    family = "gaussian", data=data,
    iterations = its, warmup = burnin, burnin = burnin,
    hyperparams=list(a=a, b=b), nuts_settings=list(adapt_delta = 0.9, max_treedepth=16), seed = seed)
  
  
  exp_ident <-res$model$exp_ident
  
  # beta samples
  betas <- res$samples$beta[(burnin+1):its, ]
  
  #beta_tilde samples
  bts <- betas
  bts[,exp_ident] <- exp(bts[,exp_ident])
  
  bt <- colMeans(bts)
  
  
  # construct grid for prediction check
  pred_grid <- dgp(pnon = pnon)


  ngrid <- 100
  ys <- seq(from = min(pred_grid$y), to = max(pred_grid$y), length = ngrid)

  # get true distribution at grid points
  pt <- pdf(truth, pred_grid, q=ys)
  
  # construct prediction basis matrices etc.
  Bpreds <- Bppreds <- vector("list", pnon+2)
  
  bctm_est <- bctm_est_med <- bctm_est_mean <- matrix(0, nrow=nrow(pred_grid), ncol=ngrid)
  X_mats <- vector("list", nrow(pred_grid))
  
  pred_grid["ng"] <- NULL

  
  for(i in 1:nrow(pred_grid)){
    pred_i <- pred_grid[i,]
    
    xsmatb<- model.matrix(mod,  data=data.frame(y=ys ,pred_i[varn])) %*%t(bts)
    bctm_est_med[i,] <- apply(pnorm(xsmatb), 1, median)
    bctm_est_mean[i,] <- apply(pnorm(xsmatb), 1, mean)
    
    bctm_est[i,] <- pnorm(model.matrix(mod,  data=data.frame(y=ys ,pred_i[varn])) %*%bt )
  }
  
  
  
  pbctms <- (colMeans(abs(pt - t(bctm_est))))
  pbctms_med <- (colMeans(abs(pt - t(bctm_est_med))))
  pbctms_mean <- (colMeans(abs(pt - t(bctm_est_mean))))
  
  
  
  # 1) true 90% quantiles for grid
  alpha_grid <- seq(0.1, 0.9, by=0.1)
  
  
  compute_alpha_deviation <- function(alpha){
    
    tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    # 2) calculate F(hhat(q_true,x)) 
    test_grid <- pred_grid
    test_grid$y <- tru_q
    
    
    
    Xpred <- model.matrix(mod, test_grid)
    
    # 3) deviation from alpha
    h_est <- Xpred%*%bt
    pnorm(h_est) - alpha
    
  }
  
  
  alpha_df_bctm <- sapply(alpha_grid, function(x) compute_alpha_deviation(x)) %>% as_tibble %>% set_names(alpha_grid)
  
  
  compute_quantile_deviation <- function(alpha){
    
    h <- function(q, data){
      
      pnon <- ncol(data)-2
      pred_data <- data.frame(y=q, data)
      Xpred <- model.matrix(mod, pred_data)
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
    
    
    q_preds - tru_q
  }
  
  quantile_df_bctm <- sapply(alpha_grid, function(x) compute_quantile_deviation(x)) %>% as_tibble %>% set_names(alpha_grid)
  
  
  
  
  tru_pq1 <- pnorm(qnorm(0.25, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)))
  tru_pq3 <-  pnorm(qnorm(0.75, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)))
  
  
  compute_quantile_ranges <- function(alpha){
    
    h <- function(q, data){
      
      pnon <- ncol(data)-2
      pred_data <- data.frame(y=q, data)
      
      Xpred <- model.matrix(mod, pred_data)
      
      h_est <- Xpred%*%bt
      h_est
    } 
    
    inverse <- function(f, data,lower, upper){
      function(y){
        uniroot(function(x){f(x, data) - y}, lower = lower, upper = upper, tol=1e-3)[1]
      }
    }
    
    q_preds <- rep(0, nrow(pred_grid))
    
    for(i in 1:nrow(pred_grid)){
      h_inverse <- inverse(h, data=pred_grid[i,1:(pnon+2)],-100, 100)
      q_preds[i] <- h_inverse(qnorm(alpha))$root
    }
    
    
    q_preds
  }
  
  
  
  alpha_borders <- c(0.25, 0.75)
  quantile_ranges_df_bctm <- sapply(alpha_borders, function(x) compute_quantile_ranges(x)) %>% as_tibble %>%  set_names(alpha_borders)
  
  
  upper <- pull(quantile_ranges_df_bctm, 2)
  lower <- pull(quantile_ranges_df_bctm, 1)
  
  
  interval_score_bctm <- interval_score(true_values = pred_grid$y,
                                        lower = lower,
                                        upper = upper,
                                        interval_range = 0.5)
  
  
  tru_q1 <- qnorm(tru_pq1)
  tru_q3 <- qnorm(tru_pq3)
  
  interval_score_tru <- interval_score(true_values = pred_grid$y,
                                       lower = tru_q1,
                                       upper = tru_q3,
                                       interval_range = 0.5)
  
  
  
  
  
  
  
  
  # ----------------------------# ----------------------------# --------------------------------------------------------------------------------------------------
  # MLT
  # --------------------------------------------------------------------------------------------------
  
  var_y <- numeric_var("y", support=range(y), bounds=range(y))
  
  B_y <- Bernstein_basis(var_y, order = 1, ui = "increasing")
  
  fm_yx <- as.formula(paste0("y ~ ", paste0(varn, collapse = " + ")))
  ctm_yx <- ctm(B_y, interacting = fm_yx[-2L], data = data,todistr = "Normal", sumconstr=F)
  mlt_yx <- mlt(ctm_yx, data = data)
  pmlt <- predict(mlt_yx, newdata=pred_grid, q=ys, type = "distribution")
  pmlts <-(colMeans(abs(pt - (pmlt))))
  summary(pmlts)  
  
  
  
  
  alpha_df_mlt <- sapply(alpha_grid, function(alpha){
    
    
    tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    sapply(1:nrow(tru_q), function(i)
      predict(mlt_yx, newdata=pred_grid[i,], q=tru_q[i], type = "distribution")  - alpha
    )
  }
  )%>% as_tibble %>% set_names(alpha_grid)
  
  
  quantile_df_mlt <- sapply(alpha_grid, function(alpha){
    
    
    tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    sapply(1:nrow(tru_q), function(i)
      predict(mlt_yx, newdata=pred_grid[i,], prob=alpha, type = "quantile")  - tru_q[i]
    )
  }
  )%>% as_tibble %>% set_names(alpha_grid)
  
  
  quantile_ranges_mlt <- sapply(alpha_borders, function(alpha){
    tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    
    sapply(1:nrow(tru_q), function(i)
      predict(mlt_yx, newdata=pred_grid[i,], prob=alpha, type = "quantile")  
    )
  }  )%>% as_tibble %>% set_names(alpha_borders)
  
  
  upper <- pull(quantile_ranges_mlt, 2)
  lower <- pull(quantile_ranges_mlt, 1)
  
  
  interval_score_mlt <- interval_score(true_values = pred_grid$y,
                                       lower = lower,
                                       upper = upper %>% unlist,
                                       interval_range = 0.5)
  
  
  
  # predict median
  ys <- seq(from = min(xdf$y), to = max(xdf$y), length = 100)
  
  pdf <- function(object, newdata, q)    UseMethod("pdf")   
  
  pdf.truth <- function(object, newdata, q = ys) {
    m <- truth$mean(newdata)
    s <- truth$sd(newdata)
    sapply(1:length(m), function(i) 
      pnorm(q, mean = m[i], sd = s[i])
    )
  }
  
  qf.truth <- function(object, newdata, p) {
    m <- truth$mean(newdata)
    s <- truth$sd(newdata)
    sapply(1:length(m), function(i) 
      qnorm(p, mean = m[i], sd = s[i])
    )
  }
  
  
 
  mad_summaries <- rbind(c(model="bctm", summary(pbctms)),  #bctm_sum_mean = summary(pbctms_mean),  bctm_sum_med = summary(pbctms_med),
                         c(model="mlt", summary(pmlts))
                        )
  
  alpha_summaries = rbind(data.frame(model="bctm", alpha_df_bctm), 
                          data.frame(model="mlt", alpha_df_mlt))
  
  quantile_summaries <- rbind(data.frame(model="bctm", quantile_df_bctm),
                              data.frame(model="mlt", quantile_df_mlt))
  
  interval_summaries <- rbind(data.frame(model="bctm" , is=(interval_score_bctm)), 
                                   data.frame(model="mlt" , is=(interval_score_mlt)))
  
  
  results <- list(MAD_summaries = mad_summaries,
                  alpha_summaries = alpha_summaries,
                  quantile_summaries = quantile_summaries, 
                  interval_summaries = interval_summaries,
                  seed = seed)
  
  return(results)
  
} 

# length of markov chain
its <- 4000
# lenght of nuts warmup
warmup <- 2000
burnin <- 2000
# number of parallel chains
nchain <- 1
# no of cores used
cores <- 1
# no  of replications
R <- 100

m <- 8

if(nchain*cores > 50) stop("Too many cores selected for parallel computation")


# no of noise parameters
pmin <- 0
pmax <- 5

#MADS <- vector("list", pmax+1)

a <- 1
b <- 0.001

# 5 x 20 cores
ret <- vector(mode = "list", length = pmax + 1)
names(ret) <- ps <- pmin:pmax

for(pnon in pmin:pmax){
  seeds <- matrix(round(runif(R) * 1000), nrow = 1)
  
  for(j in 1:nrow(seeds)){
    print(paste0("Current Scenario: pnon= ", pnon))
    cat("j: ", j, " pnon: ", pnon, "\n")
    results <- mclapply(seeds[j,], run_sims, its=its, burnin=burnin,  pnon=pnon, mc.cores = cores, m = m, a = a, b = b )          
    ret[[as.character(pnon)]]  <- c(ret[[as.character(pnon)]],
                                    results)
    
  }
    res <- ret[as.character(pnon)]
 save(res,
      file =  paste0("/processed_data/lin_sims_nonlin_pnon", pnon, "_its", its,"_burnin",burnin,"_m", m, ".RData" ))  
}

