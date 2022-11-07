##
## Script name: generate_sims_nonlin.R
##
## Purpose of script: implements nonlinear simulation scenario as shown in the paper
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------

library(dplyr)

packages <- c("tidyverse", "stringr","MASS","doParallel", "scoringutils", "BayesX", "Rcpp", "RcppArmadillo", "splines", "mgcv", "Matrix", "MCMCpack", "sdPrior", "R2BayesX", "RhpcBLASctl", "scam", "bamlss", "mlt")

source("code/helpers.R")

load_inst_packages(packages)

# library("RhpcBLASctl")
# omp_set_num_threads(1)
# blas_set_num_threads(1)


# code provided by hothorn et al. (2017)
source("code/mlt_sims/simfuns.R")


sourceCpp("code/rcpp/posterior_grad_xx2.cpp")



source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts_omega.R")
source("code/nuts/adnuts_helper.R")


run_sims <- function(seed, its, burnin,  pnon, m, a = 1, b = 0.001){
  
  set.seed(seed)
  
  data <-  dgp(c(100, 100), pnon = pnon)
  n <- nrow(data)
  
  y<- data$y #<- as.numeric(scale(df$y))
  
  varn <- colnames(data)
  varn <- varn[grep("^x", varn)]
  
  # select matrix of covariates
  xmat <- xlin <- data %>% dplyr::select(all_of(varn)) %>% as.matrix()
  ymat <- as.matrix( y)
  
  f <- paste("y", "~", paste("hyx_sm(y," , varn, ",data=data, q=c(10,10), add_to_diag=10e-6)",
                             collapse=" + "))
  
  
  res <- bctm(      
    as.formula(f),
    family = "gaussian", data=data,
    iterations = its, warmup = burnin, burnin = burnin,
    hyperparams=list(a=a, b=b), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)
  
  
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
  # ys <- data$y
  # get true distribution at grid points
  
  pt <- pdf(truth, pred_grid, q=ys)
  
  # prediction grid for y
  
  # construct prediction basis matrices etc.
  Bpreds <- Bppreds <- vector("list", pnon+2)
  
  bctm_est <- bctm_est_med <- bctm_est_mean <- matrix(0, nrow=nrow(pred_grid), ncol=ngrid)
  X_mats <- vector("list", nrow(pred_grid))
  
  # want ~ 1000 samples of estimated cdfs
  Bpreds_list <- Bppreds_list <- vector("list", pnon+2)
  
  
  for(i in 1:nrow(pred_grid)){
    xsmat <- as.matrix(pred_grid[i, colnames(pred_grid) %in% varn], ngrid, 1)
    
    # data <- expand.grid(ys, xsmat)
    
    for(j in 1:(pnon+2)){
      pred_data <- expand.grid(ys, xsmat[j])
      colnames(pred_data) <- c("y", varn[j])
      
      Bpreds_list[[j]] <- res$predvars[[j]]$B(pred_data)
      Bppreds_list[[j]] <- res$predvars[[j]]$Bp(pred_data)
      
      
    }
    
    
    Xpred <- cbind(1,do.call(cbind, Bpreds_list))
    
    X_mats[[i]] <- Xpred
    bctm_est_med[i,] <- apply(pnorm(Xpred%*%t(bts)), 1, median)
    bctm_est_mean[i,] <- apply(pnorm(Xpred%*%t(bts)), 1, mean)
    
    bctm_est[i,] <- pnorm(Xpred%*%bt )
    
    
  }
  
  
  
  pbctms <- (colMeans(abs(pt - t(bctm_est))))
  pbctms_med <- (colMeans(abs(pt - t(bctm_est_med))))
  pbctms_mean <- (colMeans(abs(pt - t(bctm_est_mean))))
  
  
  
  # 1) true 90% quantiles for grid
  alpha_grid <- seq(0.1, 0.9, by=0.1)
  
  
  compute_alpha_deviation <- function(alpha){
    
    tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    # 2) F(hhat(q_true,x)) berechnen
    test_grid <- pred_grid
    test_grid$y <- tru_q
    
    for(j in 1:(pnon+2)){
      pred_data <- test_grid[, c(paste0("x", j), "y")]
      Bpreds_list[[j]] <- res$predvars[[j]]$B(pred_data)
      
    }
    
    Xpred <- cbind(1,do.call(cbind, Bpreds_list))
    
    # 3) Abweichung von alpha bestimmen
    h_est <- Xpred%*%bt
    pnorm(h_est) - alpha
    
  }
  
  
  alpha_df_bctm <- sapply(alpha_grid, function(x) compute_alpha_deviation(x)) %>% as_tibble %>% set_names(alpha_grid)
  
  
  
  
  compute_quantile_deviation <- function(alpha){
    
    h <- function(q, data){
      
      pnon <- ncol(data)-2
      pred_data <- data.frame(y=q, data)
      for(j in 1:(pnon+2)){
        # pred_data <- test_grid[, c(paste0("x", j), "y")]
        Bpreds_list[[j]] <- res$predvars[[j]]$B(pred_data)
        
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
    
    
    q_preds - tru_q
  }
  
  

  quantile_df_bctm <- sapply(alpha_grid, function(x) compute_quantile_deviation(x)) %>% as_tibble %>% set_names(alpha_grid)
  
  

  
  tru_pq1 <- pnorm(qnorm(0.25, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)))
  tru_pq3 <-  pnorm(qnorm(0.75, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)))
  
  
  compute_quantile_ranges <- function(alpha){
    
    h <- function(q, data){
      
      pnon <- ncol(data)-2
      pred_data <- data.frame(y=q, data)
      
      for(j in 1:(pnon+2)){
        # pred_data <- test_grid[, c(paste0("x", j), "y")]
        Bpreds_list[[j]] <- res$predvars[[j]]$B(pred_data)
        
      }
      
      Xpred <- cbind(1,do.call(cbind, Bpreds_list))
      
      # 3) Abweichung von alpha bestimmen
      h_est <- Xpred%*%bt
      h_est
    } 
    
    #tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    
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
    
    
    q_preds #- tru_q
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
  
  
  compute_quantile <- function(alpha){
    
    h <- function(q, data){
      
      pnon <- ncol(data)-2
      pred_data <- data.frame(y=q, data)
      for(j in 1:(pnon+2)){
        # pred_data <- test_grid[, c(paste0("x", j), "y")]
        Bpreds_list[[j]] <- res$predvars[[j]]$B(pred_data)
        
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
  
  quantile_score_df_bctm <- sapply(c(0.05, 0.95), function(x) compute_quantile(x)) %>%
    as_tibble %>% set_names(c(0.05, 0.95))
  
  # quantile score
  calculate_quantile_score <- function(y, q, alpha){
    2*((y < q) - alpha)*(alpha - y)
    
  }
  
 bctm_quantile_scores <- c(calculate_quantile_score(pred_grid$y, quantile_score_df_bctm$`0.05`, 0.05) %>% mean,
  calculate_quantile_score(pred_grid$y, quantile_score_df_bctm$`0.95`, 0.95) %>% mean)
  
  
  
  
  
  # ----------------------------# ----------------------------# --------------------------------------------------------------------------------------------------
  # MLT
  # --------------------------------------------------------------------------------------------------
  
  var_y <- numeric_var("y", support=range(y), bounds=range(y))
  B_y <- Bernstein_basis(var_y, order = m+2, ui = "increasing")
  
  
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
  
  pmlt <- predict(mlt_yx, newdata=pred_grid, q=ys, type = "distribution")
  pmlts <-(colMeans(abs(pt - (pmlt))))
  
  
  
  
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
  
  
  quantile_scores_df_mlt <- sapply(c(0.05, 0.95), function(alpha){
    # tru_q <- qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid))
    
    
    qs <- sapply(1:nrow(pred_grid), function(i)
      predict(mlt_yx, newdata=pred_grid[i,], prob=alpha, type = "quantile")  
    )
      2*((pred_grid$y < qs) - alpha)*(alpha - pred_grid$y)
    
  }  )%>% as_tibble %>% set_names(c(0.05, 0.95))
  
  mlt_quantile_scores <- colMeans(quantile_scores_df_mlt)
  
  # mlt_devs_df <-  tibble(abs(pnorm(quantile_ranges_mlt[,1, drop=T]) - tru_pq1) , abs(pnorm(quantile_ranges_mlt[,2, drop=T]) - tru_pq3)) %>% set_names(c("PQ_1", "PQ_3"))
  
  upper <- pull(quantile_ranges_mlt, 2)
  lower <- pull(quantile_ranges_mlt, 1)
  
  
  interval_score_mlt <- interval_score(true_values = pred_grid$y,
                                       lower = lower,
                                       upper = upper %>% unlist,
                                       interval_range = 0.5)
  
  
  
  # --------------------------------------------------------------------------------------------------
  # BAMLSS - distributional regression
  # --------------------------------------------------------------------------------------------------
  
  varn <- colnames(data)
  varn <- varn[grep("^x", varn)] 
  varn <- varn[!(varn %in% c("x1", "x2"))]
  mu.fm <- "y ~ s(x = x1, by = x2)"
  # mu.fm <- "y ~ te(x1, x2, k=c(8,8))" #+te(x1, x3, k=c(8,8)) + te(x1, x4, k=c(8,8))+ te(x1, x5, k=c(8,8))"
  if (length(varn) > 0) {
    mu.fm <- paste(mu.fm, "+", paste(paste("s(", varn, ")"), 
                                     collapse = "+"))
  }
  mu.fm <- as.formula(mu.fm)
  si.fm <- as.formula(sigma ~ s(x1))
  
  f <- list(mu.fm, si.fm)
  
  
  b1 <- bamlss(f, family = "gaussian", data = data)
  mu <- predict(b1, newdata = pred_grid, type="parameter")
  
  pbamlss <-sapply(1:length(mu$mu), function(i) {
    pnorm(ys, mean = mu$mu[i], sd = mu$sigma[i])
  })
  pbamlss <- colMeans(abs(pt - pbamlss))
  summary(pbamlss)
  print("BAMLSS OK")
  
  
  
  # true quantiles for prediction grid
  tru_q <-sapply(alpha_grid, function(alpha) qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)))
  
  pred_bamlss <- predict(b1, newdata=pred_grid,  type="parameter")
  
  alpha_df_bamlss <- sapply(1:9, function(i) pnorm(tru_q[,i], mean=pred_bamlss$mu, sd=pred_bamlss$sigma) - alpha_grid[i]) %>% as_tibble %>% set_names(alpha_grid)
  
  # --------------------------------------------------------------------------------------------------
  # BAMLS - quantile regression
  # --------------------------------------------------------------------------------------------------
  
  # get true quantiles
  tru_q <- sapply(alpha_grid,
                  function(alpha) qnorm(alpha, mean=truth$mean(pred_grid), sd = truth$sd(pred_grid)
                  ))
  
  f <- "y ~ s(x1) + s(x2)"
  # mu.fm <- "y ~ te(x1, x2, k=c(8,8))" #+te(x1, x3, k=c(8,8)) + te(x1, x4, k=c(8,8))+ te(x1, x5, k=c(8,8))"
  if (length(varn) > 0) {
    f <- paste(f, "+", paste(paste("s(", varn, ")"), 
                             collapse = "+")) %>% as.formula
  }
  
  ## Quantile regression.
  bamlss_models <- lapply(alpha_grid, function(alpha) bamlss(f,
                                                             data = data, optimizer = FALSE, sampler = BayesX,
                                                             family = gF("quant", prob = alpha))
  )
  
  # prediction grid for y
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
  
  
  # probailities for all prediction quantiles for all x combinations in pred_grid
  quantile_df_bamlss_qr <- sapply(1:9, function(i){
    
    q_alpha <- predict(bamlss_models[[i]], newdata=pred_grid)
    q_alpha - tru_q[,i]
  }) %>% as_tibble %>% set_names(alpha_grid)
  
  
  bamlss_qr_q1 <- bamlss(f, data = data, optimizer = FALSE, sampler = BayesX, family = gF("quant", prob = 0.25))
  bamlss_qr_q3 <- bamlss(f, data = data, optimizer = FALSE, sampler = BayesX, family = gF("quant", prob = 0.75))
  
  q1_pred_bamlss_qr <- predict(bamlss_qr_q1, newdata=pred_grid)
  q3_pred_bamlss_qr <- predict(bamlss_qr_q3, newdata=pred_grid)
  
  
  
  interval_score_bamlss_qr <- interval_score(true_values = pred_grid$y,
                                             lower = q1_pred_bamlss_qr,
                                             upper = q3_pred_bamlss_qr,
                                             interval_range = 0.5)
  
  
  bamlss_qr_q005 <- bamlss(f, data = data, optimizer = FALSE, sampler = BayesX, family = gF("quant", prob = 0.05))
  bamlss_qr_q095 <- bamlss(f, data = data, optimizer = FALSE, sampler = BayesX, family = gF("quant", prob = 0.95))
  
  q005_pred_bamlss_qr <- predict(bamlss_qr_q005, newdata=pred_grid)
  q095_pred_bamlss_qr <- predict(bamlss_qr_q095, newdata=pred_grid)
  
  
  bamlss_qr_q_scores = c((2*((pred_grid$y < q005_pred_bamlss_qr) - 0.05)*(0.05 - pred_grid$y)) %>% mean,
  (2*((pred_grid$y < q095_pred_bamlss_qr) - 0.95)*(0.95 - pred_grid$y)) %>% mean)
  
  mad_summaries <- rbind(c(model="bctm", summary(pbctms)),
                         c(model="mlt", summary(pmlts)),
                         c(model="bamlss", summary(pbamlss)))
  
  alpha_summaries = rbind(data.frame(model="bctm", alpha_df_bctm), data.frame(model="mlt", alpha_df_mlt), data.frame(model="bamlss_distreg", alpha_df_bamlss))
  
  quantile_summaries <- rbind(data.frame(model="bctm", quantile_df_bctm),
                              data.frame(model="mlt", quantile_df_mlt),
                              data.frame(model="bamlss_qr", quantile_df_bamlss_qr))
  
  interval_summaries <- rbind(data.frame(model="bctm" , is=(interval_score_bctm)), 
                                   data.frame(model="mlt" , is=(interval_score_mlt)),
                                   data.frame(model="bamlss_qr" , is=(interval_score_bamlss_qr))
  )
  
  quantile_scores <- rbind(mlt_quantile_scores, bctm_quantile_scores, bamlss_qr_q_scores)
  

  results <- list(MAD_summaries = mad_summaries,
                  alpha_summaries = alpha_summaries,
                  quantile_summaries = quantile_summaries, 
                  interval_summaries = interval_summaries,
                  quantile_scores = quantile_scores,
                  seed = seed)
  
  return(results)
  
} 

set.seed(42)

# length of markov chain
its <- 4000
# lenght of nuts warmup
warmup <- 2000
burnin <- 0.5 * its
# number of parallel chains
nchain <- 1
# no of cores used
cores <- 20
# no  of replications
R <- 100

# no. of knots in marginal spline
m <- 8

if(nchain*cores > 50) stop("Too many cores selected for parallel computation")

# max number of noise parameters
pmax <- 5


a <- 1
b <- 0.001

# 5 x 20 cores
ret <- vector(mode = "list", length = pmax + 1)
names(ret) <- ps <- 0:pmax

# iterates through scenarios with increasing number of noise features
for(pnon in 0:pmax){
  seeds <- matrix(round(runif(R) * 1000), nrow = 1)
  
  for(j in 1:nrow(seeds)){
    print(paste0("Current Scenario: pnon= ", pnon))
    cat("j: ", j, " pnon: ", pnon, "\n")
    results <- mclapply(seeds[j,], function(seed) run_sims(its=its, burnin=burnin,  pnon=pnon, m = m, a = a, b = b, seed=seed) , mc.cores = cores)          
    ret[[as.character(pnon)]]  <- c(ret[[as.character(pnon)]],
                                    results)
    
  }
  res <- ret[as.character(pnon)]
  save(res,
       file =  paste0(mainpath,"../processed_data/new_sims_nonlin_pnon", pnon, "_its", its,"_burnin",burnin,"_m", m, ".RData" ))

}
