##
## Script name: plot_coverages.R
##
## Purpose of script: plots coverages obtained from BCTM and MLT as shown in the paper and supplement
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------

library(tidyverse)

# 
R <- 100

for(n in c(100, 500)){

load(paste0("processed_data/sims/cov_sims_n",n, ".RData"))

fig_path <- "manuscript/supplement/figs/"

res <- ret[[1]]

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
  

ggsave(filename = paste0(fig_path,"/coverages_n", n, ".pdf"), plot = p, height=4, width=7, units="in")
}