##
## Script name: plot_mads.R
##
## Purpose of script: plots the median absolute deviation results for all scenarios as shown in the paper
##
## Author: BLIND
##
## Date Created: 2020-10-7
##
## Email: BLIND
##
## ---------------------------

library(tidyverse)
library(dplyr)
library(tidyr)


pnons <- 0:5
length(pnons)
gg <-  vector("list", length(pnons))
res <- vector("list", length(pnons))

lin_MADS <- nonlin_MADS <-  all_MADS <- vector("list", length(pnons))
for(pnon in pnons){
  
  load(paste0("processed_data/sims/sims_lin_pnon" , pnon, "_its2000_burnin1000.RData"))
  lin_MADS[[pnon+1]] <- bind_rows(lapply(MADS, function(x) bind_rows(x) %>% mutate(model=c("BCTM_lin", "MLT_lin", "BAMLSS")))) %>%
    mutate(across(is.table, as.numeric)) #%>%
  
  load(paste0("processed_data/sims/score_sims_nonlin_pnon", pnon, "_its4000_burnin2000_m8.RData"))
  ret <- res
  
  res_nl <- do.call(rbind, lapply(ret[[as.character(pnon)]], "[", "MAD_summaries")) 
  res_nl <- do.call(rbind, res_nl[!is.na(res_nl) ]) %>% as_tibble
  nonlin_MADS[[pnon+1]]  <- res_nl %>% dplyr::select(c("model", "Min.", "Median", "Max.")) %>% 
    mutate_at(vars("Min.", "Median", "Max."), as.numeric) %>% 
    filter(model !="bamlss")
  
  all_MADS[[pnon+1]] <- bind_rows(lin_MADS[[pnon+1]], nonlin_MADS[[pnon+1]])
  
}


gg <- bind_rows(all_MADS, .id="pnon") %>%
  pivot_longer(cols=-c(pnon, model)) %>%
  filter(name %in% c("Min.", "Median", "Max.")) %>% 
  mutate(name = str_replace(name, "Min.", "Min. MAD"),
         name = str_replace(name, "Median", "Med. MAD"),
         name = str_replace(name, "Max.", "Max. MAD")) %>% 
  mutate(source = factor(model, levels=c("BCTM_lin", "MLT_lin", "bctm", "mlt", "BAMLSS"), labels=c("Lin. BCTM", "Lin. MLT", "Full BCTM", "Full MLT", "Oracle BAMLSS")),
         pnon = paste0("p = ",as.numeric(pnon)-1))

p1 <- gg %>% ggplot(aes(x=source, y=value)) +
  geom_boxplot() +
  # facet_wrap(name ~ pnon,ncol=length(pnons), scales="free_y")+
  facet_grid(name ~ pnon, scales="free_y") +
  theme(axis.text.x = element_text(angle = 45, size=16),
        axis.text.y = element_text(angle = 45, size=16),
        strip.background = element_rect(colour="black"),
        text = element_text(size = 18)
        ) +
  xlab("") + ylab("")
p1

ggsave("manuscript/figs/mads.png", plot=p1)

