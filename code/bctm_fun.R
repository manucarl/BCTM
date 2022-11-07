##
## Script name: bctm_fun.R
##
## Purpose of script: implements BCTM
##
## Author: BLIND
##
## Date Created: 2020-10-5
##
## Email: BLIND
##
## ---------------------------

bctm <- function(formula, family = c("gaussian", "logistic", "mev"), data, hyperparams=list(a=1, b=0.001), iterations = 2000, start=NULL,
                 warmup = iterations/2, burnin=iterations/2,  cens=NULL, seed = NULL, intercept = TRUE, looic = TRUE, 
                 nuts_settings=list(adapt_delta = 0.8, metric = NULL, step_size = NULL,
                                    adapt_mass = TRUE, max_treedepth=12, init_buffer = 75, term_buffer = 50, window = 25),
                 ...){
  
  
  family <- match.arg(family)
  
  f <- switch(as.character(family), `gaussian` = posterior_gauss,
              `logistic` = posterior_logit, `mev` = posterior_mev )
  
  ll <- switch(as.character(family), `gaussian` = ll_gauss,
              `logistic` = ll_logit, `mev` = ll_mev )
  
  grad <- switch(as.character(family), `gaussian` = gradf_gauss,
                 `logistic` = gradf_logit, `mev` = gradf_mev )
  
  
  model <- bctm_setup(formula, data, intercept)
  
  
  X <- as.matrix(model$X)
  Xp <- as.matrix(model$Xp)
  XtX <- as.matrix(t(X)%*%X)
  
  p <-  ncol(X)
  


  
  exp_inds <- model$exp_inds
  exp_ident <-   (1:ncol(X))[c(rep(F, intercept), unlist(exp_inds))]
  pen_ident <- model$pen_ident
  pen_inds <- model$pen_inds
  unpen_ident <- model$unpen_ident
  
  
  # construct model precision matrix 
  Ks <- model$Ks
  Ks_f <- unlist( model$Ks[pen_inds], recursive = F, use.names =F)
  name_groups <-names(pen_ident)
  
  labels <- model$labels
  K_inds <- lapply(labels[pen_inds], function(x) which(x == name_groups))
  
  # names(Ks_f) <- names(pen_ident)
  npen <- length(Ks_f)
  ranks <- lapply(Ks_f, rankMatrix)
  dpen_ind <- lapply(Ks,  function(x) length(x) > 1)
  spen_ind <- lapply(Ks,  function(x) length(x) == 1)
  S <- matrix(0, p - length(unlist(unpen_ident))   , p - length(unlist(unpen_ident)) )
  
  
  # model ingredients
  xx <- list(response = model$y, X = X, Xp = Xp, XtX = XtX, Xpt = t(Xp), K_inds = K_inds, p = p,
             npen = length(Ks_f) , pen_ident = model$pen_ident, unpen_ident = model$unpen_ident, m_pen_ident = unlist(unique(pen_ident))-1,
             ranks = lapply(Ks_f, rankMatrix), S = S, labels = model$labels, family = family, 
             exp_ident = exp_ident-1, n_coef =  ncol(X), hyperparams = hyperparams, Ks_f = Ks_f, Ks=Ks)
  
  
  if(!is.null(cens)){
    cind <- as.logical(cens)
    X2 <- X[!cind,]
    X1 <- X[cind,]
    
    Xp2 <- Xp[!cind,]
    Xp1 <- Xp[cind,]
    
    # if(family != "logistic") stop("currently, censoring only with logistic reference distribution.")
    f <- switch(as.character(family), `gaussian` = posterior_gauss_cens,
                `logistic` = posterior_logit_cens, `mev` = posterior_mev_cens )
    ll <- switch(as.character(family), `gaussian` = ll_gauss_cens,
                 `logistic` = ll_logit_cens, `mev` = ll_mev_cens )
    grad <- switch(as.character(family), `gaussian` = gradf_gauss_cens,
                   `logistic` = gradf_logit_cens, `mev` = gradf_mev_cens )
    
    
    xx$X1 <- X1
    xx$X2 <- X2
    xx$Xp1 <- Xp1
    xx$Xp2 <- Xp2
    
  }
  
  
  
  
  if(is.null(start)){
    message(paste("Using zeros as starting values."))
    start <- rep(0,p) 
  } else{
    message(paste("Using random draws from U[-0.5, 0.5] as starting values."))
    start <- runif(p, -0.5, 0.5)
  }
  
  chain <-NUTS(n_iter = iterations, warmup=warmup, burnin=burnin, xx=xx, f = f, gr = grad, ll=ll, start=start, nuts_settings = nuts_settings, seed = seed, theta = theta, ...)
  
  
  betas <- chain$beta[burnin:iterations,]
  colnames(betas) <- colnames(X)
  
  # traplot(betas)
  
  beta <- colMeans(betas)
  bt <- beta
  bt[exp_ident] <- exp(beta[exp_ident])
  
  
  
  bt_samples <-  betas
  bt_samples[,exp_ident] <- exp(bt_samples[,exp_ident])
  
  beta_tilde <- colMeans(bt_samples)

  h <- X%*%beta_tilde
  hp <- Xp%*%beta_tilde
  
  distr = switch(family, gaussian = function(h) pnorm(h), logistic = function(h) plogis(h), mev = function(h) bctm_pmev(h))
  dens = switch(family, gaussian = function(h, hp) dnorm(h)*hp, logistic = function(h, hp) dlogis(h)*hp, mev = function(h, hp)  bctm_dmev(h)*hp)
  
  
  
  fit <- data.frame(
  h=h,
  hp=hp,
  distr=distr(h),
  dens=dens(h,hp)
  )
  
  
  
   if(is.null(cens)){
 DIC <- {
   Dbar <- -2 * mean(unlist(chain$log_liks[burnin:iterations]), na.rm = T)
     Dhat <- -2 * sum(log(dens(h, hp)))
   pD <- Dbar - Dhat
   DIC <- pD + Dbar
   round(c(DIC = DIC, pD = pD, Dbar = Dbar, Dhat = Dhat),3)
 }
 # waic(x=t(log(dnorm(X%*%t(bt_samples))*Xp%*%t(bt_samples) )))
 # loo(x=t(dnorm(X%*%t(bt_samples), log=T) + log(Xp%*%t(bt_samples))))
 # waic(x=t((dnorm(X%*%t(bt_samples))*Xp%*%t(bt_samples) )))
 # waic(x=t(log(dnorm(X%*%t(bt_samples)) )))
 # waic(x=t(log(dnorm(X%*%t(bt_samples))*Xp%*%t(bt_samples)  )))
 # waic(x=t(log(Xp%*%t(bt_samples)  )))

  # log pointwise predicted density
  llpd <- sum(log(rowMeans(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))))

  WAIC11 <- {
      pwaic1 <- 2*sum(log(rowMeans(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))) - rowMeans(log(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))))
      WAIC1 <- -2*llpd + 2*pwaic1
     round(c(WAIC1 = WAIC1, llpd = llpd, pwaic1 = pwaic1),3)

  }

  WAIC2 <- {
    pwaic2 <- sum(apply(log(dens(X%*%t(bt_samples), Xp%*%t(bt_samples))), 1, var))
    WAIC2 <- -2*llpd + 2*pwaic2
    round(c(WAIC2 = WAIC2, llpd = llpd, pwaic2 = pwaic2), 3)

  }

 # llpd <- sum(log(rowMeans(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))))

  WAIC3 <- {
    pwaic3 <- 2*sum(1/(nrow(betas)-1)*(rowSums((log(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))) -
                      rowMeans(log(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))))^2))
    WAIC3 <- -2*llpd + 2*pwaic3
    round(c(WAIC3 = WAIC3, llpd = llpd, pwaic3 = pwaic3),3)

  }

 LOOIC <- ifelse(looic, suppressWarnings(loo::loo(t(log(dens(X%*%t(bt_samples), Xp%*%t(bt_samples)))))), NULL)
  
 } else{

   ll_i = switch(family, gaussian = dnorm, logistic = ll_logit_i, mev = ll_mev_i)

  DIC <- NULL
  WAIC1 <- loo::waic(sapply((burnin+1):iterations, function(x) ll_i(chain$beta[x,], xx)) %>% t)


  WAIC2 <- NULL
  WAIC3 <- NULL
  LOOIC <- NULL

}
  
  object <- list(y = model$y,
                 X = X,
                 Xp = Xp,
                 family = family,
                 hyper_parameters = xx$hyperparams,
                 formula = formula,
                 beta_tilde = beta_tilde,
                 fit = fit,
                 samples = chain,
                 funs = list(dens = dens, distr = distr),
                 IC = list(DIC = DIC,
                           WAIC1 = WAIC1,
                           WAIC2 = WAIC2,
                           WAIC3 = WAIC3,
                           LOOIC = LOOIC),
                 mcmc = list(iterations = iterations,
                             burnin = burnin),
                 model = list(label_inds = unlist(model$label_inds),
                              eff_inds = model$eff_ident,
                              exp_inds = exp_inds,
                              exp_ident = exp_ident,
                              specials = model$specials),
                 data = data,
                 predvars = model$predvar,
                 posterior_means = list(beta = colMeans(betas))
  )
  
  
  object
}
