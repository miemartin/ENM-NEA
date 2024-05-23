# Load libraries
require(pacman)
pacman::p_load(raster, sf, tidyverse, usdm)

# Load data

albo <- read.csv("data/albopictus_clean.csv")
head(albo)

aeg <- read.csv("data/aegypti_clean.csv")
head(aeg)

# VIF Analysis 

pnts_albo <- as.data.frame(albo)
vif.res <- vif(x = pnts_albo[,3:ncol(pnts_albo)])
vif.step <- vifstep(x = pnts_albo[,3:ncol(pnts_albo)], th = 5)
vrs <- vif.step@results$Variables %>% as.character()

pnts_aeg <- as.data.frame(aeg)
vif.res1 <- vif(x = pnts_aeg[,3:ncol(pnts_aeg)])
vif.step1 <- vifstep(x = pnts_aeg[,3:ncol(pnts_aeg)], th = 5) 
vrs1 <- vif.step1@results$Variables %>% as.character()

# Write the final table (vars selected)

pnts_albo <- dplyr::select(pnts_albo, x, y, all_of(vrs))
pnts_albo <- as_tibble(pnts_albo)
write.csv(pnts_albo, 'data/albopictus_clean_vars.csv', row.names = FALSE)

pnts_aeg <- dplyr::select(pnts_aeg, x, y, all_of(vrs1))
pnts_aeg <- as_tibble(pnts_aeg)
write.csv(pnts_aeg, 'data/aegypti_clean_vars.csv', row.names = FALSE)


# Correlation

corr1 <- cor(albo[,3:ncol(albo)], use = "pairwise.complete.obs")
corr1
write.csv(corr1, 'corr_albo.csv', row.names = FALSE)

corr2 <- cor(aeg[,3:ncol(aeg)], use = "pairwise.complete.obs")
corr2
write.csv(corr1, 'corr_aeg.csv', row.names = FALSE)
