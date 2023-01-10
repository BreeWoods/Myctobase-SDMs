## Code from Woods et al (2023) <Woods BL, Van de Putte AP, Hindell MA, Raymond B, Saunders RA, Walters A and Trebilco R (2023) 
#                                       Species distribution models describe spatial variability in mesopelagic fish
#                                       abundance in the Southern Ocean. Front. Mar. Sci. 9:981434.
#                                       doi: 10.3389/fmars.2022.981434>

## If you use this code, please cite the above paper.

## Tune model
library(gbm)
tune_grid <- expand.grid(learning.rate = c(0.05, 0.01, 0.005),
                         tree.complexity = c(1, 2, 3, 4, 5),
                         bag.fraction = c(0.5, 0.7, 0.8))
## columns to store tuning results
tune_grid$mean.null <- tune_grid$deviance.mean <- tune_grid$deviance.se <- tune_grid$n.trees <- NA_real_

## grid search
for (i in seq_len(nrow(tune_grid))) {
  ## reproducibility
  set.seed(191727)

  ## train model
  gbm.tune <- gbm.step(data = as.data.frame(abundance_df),
                       gbm.x = which(names(abundance_df) %in% pred_vars),
                       gbm.y = which(names(abundance_df) == response_var),
                       family = "gaussian",
                       learning.rate = tune_grid$learning.rate[i],
                       tree.complexity = tune_grid$tree.complexity[i],
                       bag.fraction = tune_grid$bag.fraction[i]
                       )
    ## Get stats
    tune_grid$mean.null[i] <- gbm.tune$self.statistics$mean.null
    tune_grid$deviance.mean[i] <- gbm.tune$cv.statistics$deviance.mean
    tune_grid$deviance.se[i] <- gbm.tune$cv.statistics$deviance.se
    tune_grid$n.trees[i] <- gbm.tune$gbm.call$best.trees
}
tune_grid <- tune_grid %>% mutate(deviance.expl = (mean.null - deviance.mean) / mean.null)
tune_grid <- tune_grid %>% dplyr::filter(n.trees > 1000)
write.csv(tune_grid, "tuneBRT.csv")
