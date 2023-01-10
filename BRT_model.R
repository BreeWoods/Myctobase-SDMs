## Code from Woods et al (2023) <Woods BL, Van de Putte AP, Hindell MA, Raymond B, Saunders RA, Walters A and Trebilco R (2023) 
#                                       Species distribution models describe spatial variability in mesopelagic fish
#                                       abundance in the Southern Ocean. Front. Mar. Sci. 9:981434.
#                                       doi: 10.3389/fmars.2022.981434>

## If you use this code, please cite the above paper.

### E. antarctica

## load packages
library(dplyr)
library(dismo)
library(ggplot2)

## load data
net <- readRDS("trawl_data_sector.rds")
fish <- readRDS("grps_enviro_EXTRA.rds")
## remove target trawls
fish <- fish %>% left_join(net %>% dplyr::select("EVENTID", "haul_type"), by = "EVENTID") %>%
    dplyr::filter(!grepl("Target", haul_type))
## select columns
station_cols <- c("EVENTID", "volume_filtered", "sector", "zone", "start_latitude", "start_longitude",
                  "cruise", "duration", "start_time", "SST_oiv2", "mld0", "do_200m",
                  "poc_annual_mean", "chla0", "sol_pos", "net_type", "net_depth_max",
                  "net_depth_min", "codend_mesh", "ssha_monthly_mean", "ssh_monthly_gradient",
                  "sst_oiv2_monthly_anom", "vertical_velocity_250", "bathy0")

## calculate total E. antarctica abundance for each trawl
fish_electrona <- fish %>% group_by(across(all_of(station_cols))) %>%
    dplyr::summarize(total_count = case_when(
                         any(grepl("^Electrona antarctica", scientificName)) ~ sum(individualCount[grepl("^Electrona antarctica", scientificName)]),
                         TRUE ~ 0L)) %>%
    ungroup %>%
    ## separate date and time
    mutate(Time = format(as.POSIXct(start_time, format = "%Y:%m:%d %H:%M:%S"), "%H:%M:%S"),
           Year = as.numeric(format(as.POSIXct(start_time, format = "%Y:%m:%d %H:%M:%S"), "%Y")),
           Month = as.numeric(format(as.POSIXct(start_time,format="%Y:%m:%d %H:%M:%S"), "%m"))) %>%
    dplyr::filter(!is.na(start_latitude)) ## discard anything with missing location

## select net types
fish_electrona <- fish_electrona %>% dplyr::filter(net_type %in% c("IYGPT", "RMT25", "RMT 8", "IYGPT with MIDOC")) %>%
    mutate(codend_mesh = as.factor(codend_mesh))

## abundance/presence
fish_electrona <- fish_electrona %>% mutate(organismQuantity = total_count / volume_filtered,
                                            abundance_1000m3 = organismQuantity * 1000)
fish_pa <-  fish_electrona %>% mutate(pres_abs = as.integer(total_count > 0))

cv.boot <- 10 # number of replicates (# of times models are launched/replicated)

## Gaussian model
abundance_df <- fish_electrona %>%
    ## remove NA volume filtered and absences
    dplyr::filter(total_count > 0 & !is.na(abundance_1000m3)) %>%
    mutate(logabundance = log(abundance_1000m3))

## check overall density
ggplot(abundance_df, aes(x = logabundance)) + geom_density() + theme_bw()

## go to tuneBRT.R script to determine optimal tc, lr, bf.

#saveRDS(abundance_df, file = "eantarctica-usage-gaussian_abun.rds")

## Stores the contribution
contr_n <- matrix(NA_real_, nrow = 15,ncol =  cv.boot, dimnames = list(c("SST_oiv2", "mld0", "do_200m", "poc_annual_mean", "chla0", "sol_pos", "net_type", "net_depth_max", "net_depth_min", "codend_mesh", "ssha_monthly_mean", "ssh_monthly_gradient", "sst_oiv2_monthly_anom", "vertical_velocity_250", "bathy0"), NULL))
summary <- matrix(NA_real_, nrow = 3, ncol = cv.boot, dimnames = list(c("mean.null", "deviance.mean", "deviance.se"), NULL))
response_var <- "logabundance"
pred_vars <- c("SST_oiv2", "mld0", "do_200m", "poc_annual_mean", "chla0", "sol_pos", "net_type", "net_depth_max", "net_depth_min", "codend_mesh", "ssha_monthly_mean", "ssh_monthly_gradient", "sst_oiv2_monthly_anom", "vertical_velocity_250", "bathy0")

## build model
for (j in seq_len(cv.boot)) {
    ## fit gbm
    ## set.seed(j) ## for reproducibility/testing
    fit2 <- gbm.step(data = as.data.frame(abundance_df),
                     gbm.x = which(names(abundance_df) %in% pred_vars),
                     gbm.y = which(names(abundance_df) == response_var),
                     family = "gaussian",
                     tree.complexity = 5,
                     learning.rate = 0.01,
                     bag.fraction = 0.8)

    ## Get predictor variable contributions
    rc <- summary(fit2, plotit = FALSE) ## extract the contribution
    contr_n[match(rc$var, rownames(contr_n)), j] <- rc[, "rel.inf"]
    ## Get other stats
    summary["mean.null", j] <- fit2$self.statistics$mean.null
    summary["deviance.mean", j] <- fit2$cv.statistics$deviance.mean
    summary["deviance.se", j] <- fit2$cv.statistics$deviance.se
}

write.csv(summary, "model_summary_stats.csv")

## calculate the average contributions for each predictor
contr_mean <- apply(contr_n, 1, mean)
contr_sd <- apply(contr_n, 1, sd)
contr_total <- paste(round(contr_mean, 3), round(contr_sd, 3), sep = " ± ")
names(contr_total) <- names(contr_mean)
contr_total <- data.frame(contr_total)
colnames(contr_total) <- "Contribution (%) of environmental descriptors to the model"
barplot(sort(contr_mean, decreasing=TRUE))

## save contributions
write.csv(contr_total, "avg_contribution.csv")
