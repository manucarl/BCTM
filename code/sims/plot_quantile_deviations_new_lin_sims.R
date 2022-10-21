library(tidyverse)

setwd("D:/onedrive/bctm_jasa_revision")
pmax <- 5
ret_alpha <- ret_qq <- ret_alpha_lin <- ret_qq_lin <- vector("list", pmax)
for(pnon in  0:pmax){


  
load(paste0("processed_data/sims/score_sims_nonlin_pnon", pnon, "_its4000_burnin2000_m8.RData"))
  results <- res[[1]]
  
  ret_alpha[[pnon+1]] <-   do.call(rbind, (lapply(results, function(x) unname(x["alpha_summaries"]) )))

  ret_qq[[pnon+1]] <-   do.call(rbind, (lapply(results, function(x) unname(x["quantile_summaries"]) )))

  load(paste0("processed_data/sims/lin_sims_nonlin_pnon", pnon, "_its4000_burnin2000.RData"))
  results_lin <- res[[1]]
  ret_alpha_lin[[pnon+1]] <- do.call(rbind, (lapply(results_lin, function(x) unname(x["alpha_summaries"]) )))
  ret_qq_lin[[pnon+1]] <-   do.call(rbind, (lapply(results_lin, function(x) unname(x["quantile_summaries"]) )))
  
}
lin_alpha <- lapply(ret_alpha_lin, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" ) %>% 
  mutate(model =if_else(model == "bctm", "bctm_lin", "mlt_lin"))
nl_alpha <- lapply(ret_alpha, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" )
gg_df <- bind_rows(lin_alpha, nl_alpha) %>%
  # bind_rows(.id="pnon" ) %>%
  set_names(c("pnon", "model", paste0("alpha = ", seq(0.1, 0.9, by=0.1)))) %>% 
  pivot_longer(cols=-c(pnon, model)) %>% 
  mutate(model=factor(model, levels=c("bctm_lin", "mlt_lin", "bctm", "mlt", "bamlss_distreg"), labels=c("bctm_lin", "mlt_lin","bctm", "mlt", "bamlss_oracle")))

gg_df$pnon <- paste0("p = ", as.numeric(gg_df$pnon)-1)





p1 <- gg_df  %>% ggplot() +
  geom_boxplot(aes(x=name, y=value)) +
  geom_hline(yintercept=0,col="blue", size=1) +
  geom_point(aes(x=name, y=mean_dev), col="red", data=gg_df %>% group_by(model, name) %>% summarise(mean_dev = mean(value)) %>% ungroup)+
  facet_grid(pnon~ factor(model, levels=c("bctm_lin", "mlt_lin","bctm", "mlt", "bamlss_oracle"),labels=c("Lin. BCTM", "Lin. MLT", "Full BCTM", "Full MLT", "Oracle BAMLSS"))) +
  ylab("Alpha deviation from zero") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45),
        strip.background = element_rect(colour="black"))
p1
# ggsave("figures/alpha_devs.png", plot=p1)


lin_qq <- lapply(ret_qq_lin, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" ) %>% 
  mutate(model =if_else(model == "bctm", "bctm_lin", "mlt_lin"))
nl_qq <- lapply(ret_qq, function(x) bind_rows(x[!is.na(x),])) %>% bind_rows(.id="pnon" )

gg_df_qq <- bind_rows(lin_qq, nl_qq) %>%
  set_names(c("pnon", "model", paste0("alpha =", seq(0.1, 0.9, by=0.1)))) %>% 
  pivot_longer(cols=-c(pnon, model))


gg_df_qq$pnon <- paste0("p = ", as.numeric(gg_df_qq$pnon)-1)
gg_df_qq$model <- factor(gg_df_qq$model, levels=c("bctm_lin", "mlt_lin","bctm", "mlt", "bamlss_qr"), labels=c("Lin. BCTM", "Lin. MLT", "Full BCTM", "Full MLT", "BAMLSS QR"))

p0 <- gg_df_qq %>% filter(pnon == "p = 0" | pnon == "p = 5") %>%
  ggplot() +
  geom_boxplot(aes(x=name, y=value, colour=model)) +
  facet_grid(~pnon) +
  scale_color_viridis_d(option="B", end=0.6)+
  xlab("") +
  ylab("Quantile deviation from zero") +
  theme(axis.text.x = element_text(angle = 45),
        strip.background = element_rect(colour="black"),
        text = element_text(size = 14)) +
  geom_hline(yintercept=0, linetype=2, colour="blue")
p0

# ggsave("figures/quantile_devs_small.png", plot=p0, height=5)


p2 <- gg_df_qq %>%  ggplot() + geom_boxplot(aes(x=name, y=value)) + geom_hline(yintercept=0,col="blue", size=1)+
  geom_point(aes(x=name, y=mean_dev), col="red", data=gg_df_qq %>% group_by(model, name) %>% summarise(mean_dev = mean(value)) %>% ungroup)+
  facet_grid(pnon~ model) +
  xlab("") +
  ylab("Quantile deviation from zero") +
  theme(axis.text.x = element_text(angle = 45),
        strip.background = element_rect(colour="black"))
p2
# ggsave("figures/quantile_devs.png", plot=p2)
