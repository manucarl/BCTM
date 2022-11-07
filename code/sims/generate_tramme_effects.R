##
## Script name: generate_tramme_effects.R
##
## Purpose of script: estimates nonlinear shift effect via tramME (shown in supplement)
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------

# function generating replications----------------------------------------------
library(dplyr)
library(tidyverse)
library(tramME)
library(parallel)
library(mgcv)

# function for scaling with vector output
scale_col <- function(x, m=mean(x, na.rm=F), std=sd(x, na.rm=F)){
  (x - m) / std
}

# data -----------------------------------------------------------------------
# size of training data
n <-100
# generate 4 covariates 
xm <- matrix(runif(4 * n, -2, 2), n, 4)# %>% scale

# grid for effect predictions
x_seq <- seq(-2, 2, length=200)


# 4 different functions for nonlinear shifts
ff1 <- function(x) x
ff2 <- function(x) x + ((2*x-2)^2)/5.5
ff3 <- function(x) -x + pi*sin(pi*x)
ff4 <- function(x) .5*x + 15*(dnorm((x-.2)/.5) - dnorm(x+.4))

require(devtools)
# install_version("tramME", version = "1.0.3", repos = "http://cran.us.r-project.org")

seed <- 123

ffs <- list(ff2, ff3, ff4)



ff1s <- function(x) scale(ff1(x))
ff2s <- function(x) scale(ff2(x))
ff3s <- function(x) scale(ff3(x))
ff4s <- function(x) scale(ff4(x))


# function generating replications----------------------------------------------
eff_sim <- function(seed){
  
  set.seed(seed)
  
  # y <- rnorm(n, mean=ff1(xm[,1]) + ff2(xm[,2]) + ff3(xm[,3]) + ff4(xm[,4]))
  y <- rnorm(n, mean= ff1(xm[,1]) + ff2(xm[,2]) + ff3(xm[,3]) + ff4(xm[,4]))
  
  df <- data.frame(y=y, xm[,1:4])
  colnames(df) <- c("y", paste0("x", 1:4))
  
  
  n_coef  <- 10
  tramme_yx <- BoxCoxME(y ~ s(x1, k=n_coef) +  s(x2, k=n_coef) + s(x3, k=n_coef) + s(x4, k=n_coef),
                        data=df)
  
  effs <- smooth_terms(tramme_yx)
  
  
  result <-   lapply(effs, function(x) x[,c(1,2)] %>% set_names(c("x", "eff")) %>% rowid_to_column("name") %>% mutate(name = paste0("V", name),
                                                                                                                     eff=scale_col(eff))
  )

  
  names(result) <- paste0("f", 1:4)
  
  bind_rows(result, .id="trafo")
}

#no of replications
R <- 100

# parallel cores
cores <- 1


set.seed(123)

seeds <-round(runif(R) * 10000)

# simulation
ret <- mclapply(seeds, eff_sim, mc.cores= cores )
result <- bind_rows(ret, .id="rep") %>% as_tibble
result
tramme_gg <- result %>% select(trafo, x, rep, eff) %>% set_names(c("f", "xp", "name", "value")) %>% mutate(value=scale_col(value))

tramme_gg %>% ggplot() +
  geom_line(aes(x=xp, y=scale_col(value), group=(name))) + facet_grid(~f)

save(tramme_gg, file=paste0("processed_data/tramme_effect_sim_", n, ".RData"))

