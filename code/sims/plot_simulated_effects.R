##
## Script name: plot_simulated_effects.R
##
## Purpose of script: plots estimates nonlinear shift effects as shown in the supplement 
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------
library(tidyverse)


fig_path <- "manuscript/supplement/figs/"

# plot coverage rates for BCTM and MLT ------------------------------------------------------------
R <- 100
for (n in c(100, 500)) {

load(paste0("processed_data/sims/cov_sims_n", n, ".RData"))


res <- ret[[1]]
# x <- "mlt_cov1"
mlt <- lapply(paste0("mlt_cov", 1:4), function(x) 
  lapply(res, "[[", x) %>%
    bind_rows(.id="rep") %>% 
    group_by(x) %>% 
    summarise(value = sum(cov)/100)
) %>% 
  bind_rows(.id = "name")

bctm <- lapply(lapply(res, "[[", "covs_bctm"), function(x) cbind(x, x=seq(-2, 2, length=200))) %>% 
  bind_rows(.id = "rep") %>% 
  group_by(x) %>% 
  summarise(`1` = sum(X1)/R,
            `2` = sum(X2)/R,
            `3` = sum(X3)/R,
            `4` = sum(X4)/R) %>%
  pivot_longer(cols = - x)



p <- bind_rows(BCTM=bctm, MLT=mlt, .id = "Model")  %>% 
  ggplot() + geom_line(aes(x=x, y=value, col=Model)) +
  facet_wrap(~factor(name, labels = paste0("f[", 1:4, "]")), labeller = label_parsed) +
  geom_hline(yintercept = 0.95, linetype = "longdash", size =0.6) +
  ylab("coverage rates") + 
  theme(legend.title = element_blank(),
        legend.position = c(0.93,0.05),
        legend.key=element_blank())

p

# ggsave(filename = paste0(fig_path,"/coverages_n", n, ".pdf"), plot = p, height=4, width=7, units="in")

# plot R simulated effect estimates for the 4 effects for BCTM and MLT ---------------------------------------------------------

# 4 original functions (same as Scheipl et. al 2012)
ff1 <- function(x) x
ff2 <- function(x) x + ((2*x-2)^2)/5.5
ff3 <- function(x) -x + pi*sin(pi*x)
ff4 <- function(x) .5*x + 15*(dnorm((x-.2)/.5) - dnorm(x+.4))


ff1s <- function(x) scale(ff1(x))
ff2s <- function(x) scale(ff2(x))
ff3s <- function(x) scale(ff3(x))
ff4s <- function(x) scale(ff4(x))


# function for scaling with vector output
scale_col <- function(x){
  (x - mean(x, na.rm=F)) / sd(x, na.rm=F)
}



x_pred <- seq(-2,2,  length=200)

bctm_effs <- lapply(res, "[[", "bctm_effs")


bctm_gg <- lapply(paste0("bctm_eff", 1:4) , function(x) {
  do.call(cbind, lapply(bctm_effs, "[[", x)) %>% 
    as_tibble %>% 
    mutate(xp = x_pred) %>% 
    pivot_longer(cols=-xp)}) %>% 
  set_names(paste0("f", 1:4)) %>% 
  bind_rows(.id="f") %>%
  group_by(f) %>% 
  mutate(value=scale_col(value))



mlt_effs <- lapply(res, "[[", "mlt_effs")

mlt_gg <- lapply(paste0("mlt_eff", 1:4) , function(x) {
  do.call(cbind, lapply(lapply(mlt_effs, "[[", x), function(x) x[10,])) %>% 
    as_tibble %>% 
    mutate(xp = x_pred) %>% 
    pivot_longer(cols=-xp)}) %>% 
  set_names(paste0("f", 1:4)) %>% 
  bind_rows(.id="f") %>%
  group_by(f) %>% 
  mutate(value=-scale_col(value))


load(paste0("processed_data/sims/tramme_effect_sim_", n , ".RData"))

gg <- bind_rows(BCTM = bctm_gg, MLT = mlt_gg, tramME = tramme_gg, .id = "model")

alpha <- 0.3
p <- gg %>%
  filter(f=="f1") %>%
  ggplot() +
  geom_line(aes(x=xp, y=value, group=name), alpha=alpha)+
  facet_wrap(~model) +
  ylab(expression(f[1](x[1]))) + xlab(expression(x[1]))+
  stat_function(fun = ff1s, col ="red", size=1)
p

ggsave(filename = paste0(fig_path,"coverage_f1_n",n, "_wtramme.pdf"), plot = p, height=4, width=8, units="in")


p <- gg %>%
  filter(f=="f2") %>%
  ggplot() +
  geom_line(aes(x=xp, y=value, group=name), alpha=alpha)+
  facet_wrap(~model) +
  ylab(expression(f[2](x[2]))) + xlab(expression(x[2]))+
  stat_function(fun = ff2s, col ="red", size=1)
p

ggsave(filename = paste0(fig_path,"coverage_f2_n",n, "_wtramme.pdf"), plot = p, height=4, width=8, units="in")

p <- gg %>%
  filter(f=="f3") %>%
  ggplot() +
  geom_line(aes(x=xp, y=value, group=name), alpha=alpha)+
  facet_wrap(~model) +
  ylab(expression(f[3](x[3]))) + xlab(expression(x[3]))+
  stat_function(fun = ff3s, col ="red", size=1)
p
ggsave(filename = paste0(fig_path,"coverage_f3_n",n, "_wtramme.pdf"), plot = p, height=4, width=8, units="in")


p <- gg %>%
  filter(f=="f4") %>%
  ggplot() +
  geom_line(aes(x=xp, y=value, group=name), alpha=alpha)+
  facet_wrap(~model) +
  ylab(expression(f[4](x[4]))) + xlab(expression(x[4]))+
  stat_function(fun = ff4s, col ="red", size=1)
p
ggsave(filename = paste0(fig_path,"coverage_f4_n",n, "_wtramme.pdf"), plot = p, height=4, width=8, units="in")
}