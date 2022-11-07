##
## Script name: plot_quantile_scores.R
##
## Purpose of script: plots quantile scores as shown in the paper
##
## Author: BLIND
##
## Date Created: 2022-15-7
##
## Email: BLIND
##
## ---------------------------

library(tidyverse)
library(dplyr)
library(tidyr)



pmax <- 5
pnons <- 0:pmax
length(pnons)
gg <-  vector("list", length(pnons))
res <- vector("list", length(pnons))

quantile_score_list <- vector("list", length(pnons))
for(pnon in pnons){
  load(paste0("processed_data/sims/score_sims_nonlin_pnon", pnon, "_its4000_burnin2000_m8.RData"))
  
  quantile_scores <-  do.call(rbind, lapply(res[[as.character(pnon)]], "[", "quantile_scores"))
  
  
  quantile_scores <- do.call(rbind, quantile_scores[!is.na(quantile_scores)]) %>% as.data.frame %>%  rownames_to_column("model") 
  quantile_scores$model <- str_extract(quantile_scores$model, "[^_]+")
  
  
  quantile_score_list[[pnon+1]] <- quantile_scores
  
}

p1 <- bind_rows(quantile_score_list, .id="pnon") %>%
  mutate(model =  factor(model, levels=c("bctm", "mlt", "bamlss"), labels=c("Full BCTM", "Full MLT", "BAMLSS QR"))) %>%
  rename(Model = model) %>% 
  pivot_longer(cols=-c(pnon, Model)) %>%
  ggplot() + geom_boxplot(aes(x=name, y=value, colour=Model)) +
  facet_grid(~ factor(pnon, labels=paste0("p = ", 0:5))) +
  scale_color_viridis_d(option="B", end=0.6)+
  ylab("Quantile score") + xlab("Quantile") +
  theme(strip.background = element_rect(colour="black"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, size=20),
        axis.text.y = element_text(angle = 45, size=20))

p1

ggsave("manuscript/figs/quantile_scores.png", plot=p1, height=6)
