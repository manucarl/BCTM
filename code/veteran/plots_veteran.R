## ---------------------------
##
## Script name:  plots_veteran.R
##
## Purpose of script: creates VA density plots with added treatment variable in "Bayesian Conditional Transformation Models"
##
## Author: BLIND
##
## Date Created: 2020-10-1
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

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "caret",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

sourceCpp("code/rcpp/posterior_grad_xx2.cpp")


seed <- 123
data(cancer, package="survival")

data2 <- veteran[!as.logical(veteran$prior),]

dummies <- predict(dummyVars(~ celltype, data = data2), newdata = data2)
colnames(dummies) <- c("squam", "small", "adeno", "large")
data2 <- cbind(data2[,c("time", "status", "karno", "trt")], dummies)

# can scale to unit interval for faster sampling
data2$trt <- ifelse(data2$trt == 2, 1, 0)

family <- "logistic"


its <- 4000
n <- nrow(data)

# model 1 - po with censoring-------------------------------------------------------------------------------------------------------------------------
object <- bctm(time ~ hy_sm(time,data=data2, q=22) + hx_lin(karno) + hx_lin(adeno) + hx_lin(small) + hx_lin(squam) + hx_lin(trt),  
               family = "logistic", data=data2,
               cens=as.logical(data2$status), iterations = its, intercept=T,
               hyperparams=list(a=1, b=0.001),
               nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)


burnin <- object$mcmc$burnin


# indices of exponentiated coefficients
exp_ident <- object$model$exp_ident

# beta samples
beta_samples <- object$samples$beta[(burnin+1):its,]

# beta_tilde samples
bt_samples <- beta_samples
bt_samples[,exp_ident] <- exp(beta_samples[,exp_ident])


# posterior mean beta_tilde
bt <- colMeans(bt_samples)


# prediction time grid
time_grid <- seq(0, 500, length= 100)

# prediction design matrix (and derivative) constructors
pred_B <- lapply(object$predvars, "[[", "B")
pred_Bp <- lapply(object$predvars, "[[", "Bp")

B0pred <- pred_B[[1]](list(time=time_grid[-1]))
Bp0pred <- pred_Bp[[1]](list(time=time_grid[-1]))

bts <-lapply(object$model$eff_inds, function(x) bt[x])


# Plotting -----------------------------------------------------------------------------------------------------------------------------
object$model$label_inds

bts <-lapply(object$model$eff_inds[-1], function(x) bt[x])


# estimated transformation function
hy <- B0pred%*%bts[["hy_sm(time)"]]

# derivative of estimated transformation funtion
hyp <- Bp0pred%*%bts[["hy_sm(time)"]]

# karno_pred <- rescale(c(40,60,80), to = c(0,1), from= range(veteran$karno))
karno_pred <- c(40,60,80)

# Panel A -----------------------------------------------------------------------------------------------------------------------------

panelA_dat_notrt <- tibble(time=time_grid, 
                           "PS=40"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"] + karno_pred[1]*bt["hx_lin(karno)"])*hyp),
                           "PS=60"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"] + karno_pred[2]*bt["hx_lin(karno)"])*hyp),
                           "PS=80"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"]+ karno_pred[3]*bt["hx_lin(karno)"])*hyp))

panelA_dat_trt <- tibble(time=time_grid, 
                         "PS=40"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"] + karno_pred[1]*bt["hx_lin(karno)"] + bt["hx_lin(trt)"])*hyp),
                         "PS=60"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"] + karno_pred[2]*bt["hx_lin(karno)"]+ bt["hx_lin(trt)"])*hyp),
                         "PS=80"= c(0,dlogis(bt["int"]+hy + bt["hx_lin(small)"]+ karno_pred[3]*bt["hx_lin(karno)"]+ bt["hx_lin(trt)"])*hyp))


panelA_dat <- bind_rows(notrt = panelA_dat_notrt, trt = panelA_dat_trt, .id="trt") %>% 
  pivot_longer(cols=-c(time, trt))


pA <- panelA_dat %>% ggplot() +
  geom_line(aes(x=time, y=value, col=name, linetype=trt ), size=0.7)+
  theme_bw()+
  theme(legend.justification = c(1, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.9, 0.9)) +
  guides(linetype = FALSE)+
  labs(colour="Score")+
  xlab("Time")+ylab("Density")

# pA


# Panel B -----------------------------------------------------------------------------------------------------------------------------
panelB_dat_notrt <- tibble(time=time_grid, 
                           "large"= c(0,dlogis(bt["int"]+ hy + karno_pred[2]*bt["hx_lin(karno)"])*hyp),
                           "adeno"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(adeno)"] + karno_pred[2]*bt["hx_lin(karno)"])*hyp),
                           "small"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(small)"] + karno_pred[2]*bt["hx_lin(karno)"])*hyp),
                           "squam"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(squam)"]+ karno_pred[2]*bt["hx_lin(karno)"])*hyp))   
panelB_dat_trt <- tibble(time=time_grid, 
                         "large"= c(0,dlogis(bt["int"]+hy + karno_pred[2]*bt["hx_lin(karno)"]+ bt["hx_lin(trt)"])*hyp),
                         "adeno"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(adeno)"] + karno_pred[2]*bt["hx_lin(karno)"] + bt["hx_lin(trt)"])*hyp),
                         "small"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(small)"] + karno_pred[2]*bt["hx_lin(karno)"]+ bt["hx_lin(trt)"])*hyp),
                         "squam"= c(0,dlogis(bt["int"]+hy+ bt["hx_lin(squam)"]+ karno_pred[2]*bt["hx_lin(karno)"]+ bt["hx_lin(trt)"])*hyp))  

panelB_dat <- bind_rows(notrt = panelB_dat_notrt, trt = panelB_dat_trt, .id="trt") %>% 
  pivot_longer(cols=-c(time, trt))

pB <- panelB_dat %>% ggplot() +
  geom_line(aes(x=time, y=value, col=name, linetype=trt ), size=0.7)+
  xlab("Time (days)") +
  ylab("Density")+
  guides(color=guide_legend(title="Cancer cell type"),
         linetype = FALSE)+
  theme_bw()+
  theme(legend.justification = c(1, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.9, 0.9))

# pB

p <- cowplot::plot_grid(pA, pB, labels=c("A","B"))
p
ggsave("manuscript/figs/vet_densities.pdf", plot =p, height=4, width=8, units="in", bg="transparent")


