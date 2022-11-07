##
## Script name: plot_quantile_deviations.R
##
## Purpose of script: plots quantile deviations as shown in the papers
##
## Author: BLIND
##
## Date Created: 2022-12-7
##
## Email: BLIND
##
## ---------------------------

library(tidyverse)

pmax <- 5
ret_alpha <- ret_qq <- vector("list", pmax)
for(pnon in  0:pmax){
  load(paste0("processed_data/sims/score_sims_nonlin_pnon", pnon, "_its4000_burnin2000_m8.RData"))
  results <- res[[1]]
  
  ret_alpha[[pnon+1]] <-   do.call(rbind, (lapply(results, function(x) unname(x["alpha_summaries"]) )))
  
  ret_qq[[pnon+1]] <-   do.call(rbind, (lapply(results, function(x) unname(x["quantile_summaries"]) )))
  
}


gg_df <- lapply(ret_alpha, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" ) %>%
  set_names(c("pnon", "model", paste0("alpha = ", seq(0.1, 0.9, by=0.1)))) %>% 
  pivot_longer(cols=-c(pnon, model)) %>% 
  mutate(model=factor(model, levels=c("bctm", "mlt", "bamlss_distreg"), labels=c("bctm", "mlt", "bamlss_oracle")))

gg_df$pnon <- paste0("p = ", as.numeric(gg_df$pnon)-1)


# alpha deviation plots
p1 <- gg_df  %>% ggplot() +
  geom_boxplot(aes(x=name, y=value)) +
  geom_hline(yintercept=0,col="blue", size=1) +
  geom_point(aes(x=name, y=mean_dev), col="red", data=gg_df %>% group_by(model, name) %>% summarise(mean_dev = mean(value)) %>% ungroup)+
  facet_grid(pnon~ factor(model, labels=c("Full BCTM", "Full MLT", "Oracle BAMLSS"))) +
  ylab("Alpha deviation from zero") +
  xlab("") +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, size=20),
        axis.text.y = element_text(angle = 45, size=20),
        strip.background = element_rect(colour="black"))
p1
ggsave("manuscript/figs/alpha_devs.png", plot=p1)


# quantile deviation plot for paper
gg_df_qq <- lapply(ret_qq, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" ) %>%
  set_names(c("pnon", "model", paste0("alpha =", seq(0.1, 0.9, by=0.1)))) %>% 
  pivot_longer(cols=-c(pnon, model))


gg_df_qq$pnon <- paste0("p = ", as.numeric(gg_df_qq$pnon)-1)
gg_df_qq$model <- factor(gg_df_qq$model, labels=c("Full BCTM", "Full MLT", "BAMLSS QR"))

p0 <- gg_df_qq %>% filter(pnon == "p = 0" | pnon == "p = 5") %>% mutate(Model = model) %>% 
  ggplot() +
  geom_boxplot(aes(x=name, y=value, colour=Model)) +
  facet_grid(~pnon) +
  scale_color_viridis_d(option="B", end=0.6)+
  xlab("") +
  ylab("Quantile deviation from zero") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, size=20),
        axis.text.y = element_text(angle = 45, size=20),
        strip.background = element_rect(colour="black")) +
  geom_hline(yintercept=0, linetype=2, colour="blue")
p0

ggsave("manuscript/figs/quantile_devs_small.png", plot=p0, height=5)


# quantile deviation plots for all scenarios
p2 <- gg_df_qq %>%  ggplot() + geom_boxplot(aes(x=name, y=value)) + geom_hline(yintercept=0,col="blue", size=1)+
  geom_point(aes(x=name, y=mean_dev), col="red", data=gg_df_qq %>% group_by(model, name) %>% summarise(mean_dev = mean(value)) %>% ungroup)+
  facet_grid(pnon~ factor(model, labels=c("Full BCTM", "Full MLT", "BAMLSS QR"))) +
  xlab("") +
  ylab("Quantile deviation from zero") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, size=20),
        axis.text.y = element_text(angle = 45, size=20),
        strip.background = element_rect(colour="black"))
p2
ggsave("manuscript/figs/quantile_devs.png", plot=p2)
