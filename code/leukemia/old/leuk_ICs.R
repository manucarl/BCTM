library(dplyr)


load("processed_data/leukemia/leuk_ph.RData")
object_ph_cens$IC$WAIC1


load("processed_data/leukemia/leuk_ph_re.RData")
object_ph_re$IC$WAIC1$estimates

load("processed_data/leukemia/leuk_ph_spat.RData")
object_ph_spat$IC$WAIC1$estimates


load("processed_data/leukemia/leuk_nph.RData")
object_nph$IC$WAIC1$estimates

mcmcplots::mcmcplot(object_nph$samples$beta)
