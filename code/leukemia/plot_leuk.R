## ---------------------------
##
## Script name: plot_leuk.R
##
## Purpose of script: generates survivor functions and spatial plot
##
## Author: BLIND
##
## Date Created: 2022-6-5
##
## Email: BLIND
##
## ---------------------------
library(dplyr)


source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts.R")
source("code/nuts/adnuts_helper.R")

packages <- c("spBayesSurv", "survival", "BayesX", "scam", "Matrix", "Rcpp", "RcppArmadillo", "MCMCpack", "sf", "rgeos", "ggplot2","RhpcBLASctl", "loo")
load_inst(packages)
# 
# omp_set_num_threads(1)
# blas_set_num_threads(1)

sourceCpp("code/rcpp/posterior_grad_xx.cpp")

# get data
data("LeukSurv", package="spBayesSurv")
data <- LeukSurv[order(LeukSurv$district), ]
n <- nrow(data) 

# function that scales prediction variable according to the scaling in training data
scale_pred <- function(x, var){
  (x - attr(var, "scaled:center"))/
    attr(var, "scaled:scale")
}

#get boundary file
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))

# construct map in graph format for neighbour matrix
nmat <- bnd2gra(nwengland)


# scale for faster sampling
data$age <- scale(data$age)#%>% as.vector
data$wbc <- scale(data$wbc)#%>% as.vector
data$tpi <- scale(data$tpi)# %>% as.vector
data$sex[data$sex==0] <- -1

seed <- 1
# no of mcmc iterations
its <- 4000


object_ph_spat <- bctm(time ~hy_sm(time, data=data,  center=T, q = 20, add_to_diag=10e-6)+
                         hx_lin(age) + hx_lin(sex)+ hx_lin(wbc) + hx_lin(tpi) +
                         hx_spat(district, nmat=nmat, data),
                       cens = as.logical(data$cens),
                       family = "mev", data=data, iterations = its, intercept=F, 
                       hyperparams=list(a=1, b=0.001), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)
# save(object_ph_spat, file="processed_data/leukemia/leuk_ph_spat_cens_plot.RData")

load("processed_data/leukemia/leuk_ph_spat_cens_plot.RData")

object <- object_ph_spat

# get model design matrix
X <- object$X
Xp <- object$Xp

burnin <- object$mcmc$burnin

# beta samples
betas <- object$samples$beta[burnin:its,]
colnames(betas) <- colnames(X)


# indicators for exponentiation
exp_ident <- object$model$exp_ident


# beta_tilde samples
bt_samples <-  betas
bt_samples[,exp_ident] <- exp(bt_samples[,exp_ident])
bt <- object$beta_tilde

# posterior means
beta_tilde <- colMeans(bt_samples)

# map regions
nwenglandsp <- bnd2sp(nwengland)
newenglandsf <-st_as_sf(nwenglandsp)

label_inds <- object$model$label_inds



# create plots -----------------------------------------------------------------------------------------------------------------
time_grid <- seq(0, 3000, length=300)

# adapt pred variables that are standardized in training data 
age_pred <-scale_pred(49, data$age)
wbc_pred <- scale_pred(38.6, data$wbc)

# construct design matrix for baseline transformation
pred_B <- lapply(object$predvars, "[[", "B")
B0pred <- pred_B[[1]](list(time=time_grid))


# get samples of transformation function for each quantile = 0.05, 0.5, 0.95 of tpi
effsamples095 <- cbind(
  h0=B0pred,
  age = age_pred,
  sex = -1,
  wbc = wbc_pred,
  tpi = quantile(data$tpi, 0.95)
) %*% t(bt_samples[, which(label_inds != "hs(district)")]) 

effsamples050 <- 
  cbind(
    B0pred,
    age = age_pred,
    sex = -1,
    wbc = wbc_pred,
    tpi = quantile(data$tpi, 0.5)
  ) %*% t(bt_samples[, which(label_inds != "hs(district)")])


effsamples005 <-
  cbind(
    h0=B0pred,
    age = age_pred,
    sex = -1,
    wbc = wbc_pred,
    tpi = quantile(data$tpi, 0.05)
  ) %*% t(bt_samples[, which(label_inds != "hs(district)")]) 

effsamples005[1,] <- effsamples050[1,] <- effsamples095[1,] <- -Inf


# minimum-extreme value distribution function
pmev <- function(x) 1 - exp(-exp(x))

# get data for plot A

gg <-
  lapply(list(
    "5%" = effsamples005,
    "50%" = effsamples050,
    "95%" = effsamples095
  ), function(x) {
    tibble(
      "pq975" = 1 - pmev(apply(x, 1, quantile, 0.975)),
      "pqmed" = 1 - pmev(apply(x, 1,  median)),
      "pq025" = 1 - pmev(apply(x, 1, quantile, 0.025)),
      time_grid
    )
  }) %>%
  bind_rows(.id = "q")

# panel A: survivor functions
p1 <- gg %>%
  ggplot(aes(x = time_grid, y = pqmed, col = q)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = pq025, ymax = pq975),
              alpha = 0.1,
              linetype = 0) +
  ylab("Survival") + xlab("Time") +
  theme(
    legend.justification = c(1, 1),
    legend.position = c(0.95, 0.95),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent")
  ) +
  scale_color_viridis_d(option = "C", begin=0.9, end=0.1) +
  labs(color = "Townsend score (quantiles)")

p1


# panel B: map ----------------------------------------------
library(viridis)
p2 <- ggplot(data = newenglandsf) +
  theme_void() +
  geom_sf(aes(fill = c(bt[label_inds == "hs(district)"]))) +
  scale_fill_viridis(option = "A",
                     begin = 1,
                     end = 0.3) +
  theme(
    legend.justification = c(1, 1),
    legend.position = c(0.99, 0.99),
    panel.border = element_rect(colour = "white", fill = NA),
    legend.box.background = element_rect(color = NA, colour = "white"),
    legend.background = element_rect(color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"),
    rect = element_rect(fill = "transparent") # all rectangles
  ) +
  labs(fill = "") #+ panel_border(remove = TRUE, color="white")
p2

p <-cowplot::plot_grid(p1, p2, labels=c("A","B"), ncol=2)#+ panel_border(remove = TRUE, color="white")
p
#ggsave("manuscript/figs/leuk_ph.pdf", plot=p,height=3.5, width=9, units="in", bg = "transparent")


