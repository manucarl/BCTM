

load_inst_packages <-function(x) {
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

lapply(packages, library, character.only = TRUE) %>%
  invisible()
}


getcoef <- function(pattern, names=NULL){
  if(is.null(names)) names <- colnames(X)
  grepl(pattern, names)
  
}



TildeFun <- function(x, ind_exp){
  res <- x
  res[ind_exp] <- exp(x[ind_exp])
  return(res)
}


# MH update of smoothing variance
tau2_update <- function(par, prop_u, u, scale_param, Ki, rank_K){
  prob <- (lp_u(u=prop_u, param=par, scale=scale_param, Ki=Ki, rank_K=rank_K) -
             lp_u(u=u, param=par, scale=scale_param, Ki=Ki, rank_K=rank_K))
  
  if (log(runif(1)) < prob){
    tau2 = exp(prop_u)
    return(tau2)
  }else{
    tau2 = exp(u)
    return(tau2)
  }
}