## Code from Woods et al (2023) <Woods BL, Van de Putte AP, Hindell MA, Raymond B, Saunders RA, Walters A and Trebilco R (2023) 
#                                       Species distribution models describe spatial variability in mesopelagic fish
#                                       abundance in the Southern Ocean. Front. Mar. Sci. 9:981434.
#                                       doi: 10.3389/fmars.2022.981434>

## If you use this code, please cite the above paper.

## Predict
library(dplyr)
library(raster)

stg_date3 <- data.frame(date = seq(as.Date("1991-01-01"), length = 348, by = "1 month"))$date
stg_date4 <- data.frame(date=seq(as.Date("1993-01-01"), length = 324, by = "1 month"))$date

## Make an exmpty raster
lims <- c(-180, 180, -75, -45)
grd <- raster(extent(lims), resolution = c(0.1, 0.1))

## read in environmental data layers
## read in mask
msk <- raster("land_ice_mask.nc")
msk <- raster("msk_land_ice.nc")
names(msk) <- "msk"
msk <- crop(msk, lims)

## sst
sst1 <- crop(stack("sea_surface_temperature_monthly_1991-2019.nc"), lims)
names(sst1) <- stg_date3

## dissolved oxygen
oxy_200 <- crop(stack("dissolved_oxygen_monthlyclimatology_Jan-Dec.nc"), lims)
names(oxy_200) <- c("do_200m_jan", "do_200m_feb", "do_200m_mar", "do_200m_apr", "do_200m_may", "do_200m_jun",
                    "do_200m_jul", "do_200m_aug", "do_200m_sep", "do_200m_oct", "do_200m_nov", "do_200m_dec")

## bathymetry
bathy0 <- raster("bathymetry.nc")
names(bathy0) <- "bathy0"
bathy0 <- crop(bathy0, grd)

## vertical velocity at 250 m
vertical_velocity_250 <- crop(raster("vertical_velocity_250.nc"), lims)

## chla
chla <- crop(stack("chla_monthly.nc"), lims)
files <- read.csv("chla_dates.csv")
names(chla) <- files$date

## mixed layer depth
mld0 <- crop(stack("mld_monthly_climatology_Jan-Dec.nc"), lims)
names(mld0) <- c("mld_200m_jan", "mld_200m_feb", "mld_200m_mar", "mld_200m_apr", "mld_200m_may", "mld_200m_jun",
                 "mld_200m_jul", "mld_200m_aug", "mld_200m_sep", "mld_200m_oct", "mld_200m_nov", "mld_200m_dec")

## poc
poc_annual_mean <- crop(raster("poc_annual_mean.nc"), lims)

## sst_anomaly
sst_oiv2_monthly_anom <- stack("sst_oiv2_monthly_anom.nc")
## rescale and mask
sst_oiv2_monthly_anom <- raster::resample(sst_oiv2_monthly_anom, grd, method = "bilinear")
sst_oiv2_monthly_anom <- raster::mask(sst_oiv2_monthly_anom, msk)

## ssh_monthly_gradient
ssh_monthly_gradient <- stack("ssh_monthly_gradient.nc")
## rescale and mask
ssh_monthly_gradient <- raster::resample(ssh_monthly_gradient, grd, method = "bilinear")
ssh_monthly_gradient <- raster::mask(ssh_monthly_gradient, msk)
names(ssh_monthly_gradient) <- stg_date4

##ssha_monthly_mean
ssha_monthly_mean <- stack("ssha_monthly_mean.nc")
## rescale and mask
ssha_monthly_mean <- raster::resample(ssha_monthly_mean, grd, method = "bilinear")
ssha_monthly_mean <- raster::mask(ssha_monthly_mean, msk)
names(ssha_monthly_mean) <- stg_date4

## prediction
mx0 <- stack(bathy0, poc_annual_mean,vertical_velocity_250)
year_preds <- lapply(1997:2011, function(y) {
    ystack <- raster::stack(lapply(c(1, 2, 11, 12), function(m) {
        doi <- paste0(if (m > 6) y - 1 else y, "-", m, "-01")
        ## monthly env data
        sstidx <- which(stg_date3 == as.Date(doi))
        mx <- stack(mx0, sst1[[sstidx]])
        names(mx)[length(names(mx))] <- "SST_oiv2"
        chlidx <- which(stg_date3 == as.Date(doi))
        mx <- stack(mx, chla[[chlidx]])
        names(mx)[length(names(mx))] <- "chla0"
        mx <- stack(mx, oxy_200[[m]])
        names(mx)[length(names(mx))] <- "do_200m"
        mx <-stack(mx, mld0[[m]])
        names(mx)[length(names(mx))] <- "mld0"
        sstaidx <- which(stg_date3 == as.Date(doi))
        mx <- stack(mx, sst_oiv2_monthly_anom[[sstaidx]])
        names(mx)[length(names(mx))] <- "sst_oiv2_monthly_anom"
        sshidx <- which(stg_date3 == as.Date(doi))
        mx <- stack(mx, ssh_monthly_gradient[[sshidx]])
        names(mx)[length(names(mx))] <- "ssh_monthly_gradient"
        sshaidx <- which(stg_date3 == as.Date(doi))
        mx <- stack(mx, ssha_monthly_mean[[sshaidx]])
        names(mx)[length(names(mx))] <- "ssha_monthly_mean"
        ## convert to data frame
        this_dx <- as.data.frame(mx, xy = TRUE)
        ## add fixed predictors (these are held fixed over the spatial grid)
        this_dx <- data.frame(this_dx, sol_pos = -30, net_depth_min = 0, net_depth_max = 200,
                              net_type = "IYGPT with MIDOC", codend_mesh = "0.5")
        saveRDS(this_dx, "gridded.data_ea_gaussian_night_IYGPT.rds")

        ## predict with model
        this_px <- predict(fit2, newdata = this_dx, type = "response")
        this_px <- exp(this_px + 1/2*sd(fit2$residuals, na.rm = FALSE)^2)
        out <- mx[[1]]
        values(out) <- this_px
        names(out) <- month.abb[m]
        out
    }))
    mean(ystack)
})

year_preds_stack <- stack(year_preds)
year_preds_stack2 <- mask(year_preds_stack, msk)
year_mean_gaus_night <- mean(year_preds_stack2)
writeRaster(year_mean_gaus_night, file = "pred_gaus_night_RMT25_0_200m.nc")

