## Re-analysis of Kronmuller & Barr (2007)

devtools::load_all("clusterperm")
devtools::load_all("exchangr")
library("tidyverse")

set.seed(1000) # for reproducibility

nmc <- 1000L # number of monte carlo runs

mod_form <- TAS ~ Speaker * Precedent * Load + Error(SubjID)

## calculate statistics on the original data
dat2 <- aov_by_bin(kb07bins, bin, mod_form)
orig <- detect_clusters_by_effect(dat2, effect, bin, stat, p)

#########################
## Effects related to Speaker
##   Speaker, Speaker:Precedent, Speaker:Load, Speaker:Precedent:Load

## need to 'fold' the data first
dat_spkr <- nest(kb07bins, -SubjID, -Speaker)

nhds_spkr <- cluster_nhds(nmc, dat_spkr, bin, mod_form,
                          shuffle_each, Speaker, SubjID)

results_spkr <- pvalues(orig, nhds_spkr) %>%
  filter(grepl("Speaker", effect))

#########################
## Effects related to Precedent
##   Precedent, Precedent:Load

dat_prec <- nest(kb07bins, -SubjID, -Precedent)

nhds_prec <- cluster_nhds(nmc, dat_prec, bin, mod_form,
                          shuffle_each, Precedent, SubjID)

results_prec <- pvalues(orig, nhds_prec) %>%
  filter(effect %in% c("Precedent", "Precedent:Load"))

#########################
## Effects related to Load
##   Load

dat_load <- nest(kb07bins, -SubjID, -Load) %>%
  arrange(SubjID)

nhds_load <- cluster_nhds(nmc, dat_load, bin, mod_form,
                          shuffle_each, Load, SubjID)

results_load <- pvalues(orig, nhds_load) %>%
  filter(effect == "Load")

#########################
## Combine all results
results <- dplyr::bind_rows(results_spkr, results_prec, results_load)
write_csv(results, "cluster_permutation_results.csv")

results
