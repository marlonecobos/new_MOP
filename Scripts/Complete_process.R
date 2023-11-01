# ------------------------------------------------------------------------------
# Project: Expanding interpretation of results from the MOP metric
# Author: Marlon E. Cobos, Hannah Owens, Jorge Soberón, and A. Townsend Peterson
# Date: 11/01/2023
# ------------------------------------------------------------------------------


# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses and figures for this project
# 
# All data required can be obtained using the code in this script. Make sure to 
# have internet connection to download such data.
#
# The example of a species accessible area used here was obtained from 
# Machado-Stredel et al. (2021; https://doi.org/10.21425/F5FBG48814) 
#
# Note: some results and figures will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages and functions needed ----------------------------------------------
# these lines load packages. If needed, install packages using
# install.packages("package_name")
library(geodata)
library(ks)
library(biosurvey)
library(mop)
library(scales)

# functions needed
source("Scripts/functions.R")
# ------------------------------------------------------------------------------


# Working directory ------------------------------------------------------------
setwd("YOUR/DIRECTORY")
# ------------------------------------------------------------------------------


# Example data -----------------------------------------------------------------
# directory for data
dir.create("Data")

# download and read accessible area simulated for the species A. ultramarina
download.file(url = "https://github.com/marlonecobos/new_MOP/raw/main/Data/M_simulated.gpkg", 
              destfile = "Data/M_simulated.gpkg", method = "wget")

au_area <- vect("Data/M_simulated.gpkg")

# downloading environmental variables 
var_current <- worldclim_global(var = "bio", res = 10, path = "Data")

## renaming variables 
names(var_current) <- gsub("wc2.1_10m_", "", names(var_current))

## selecting a subset of variables
selected <- c(5:7, 13:15)

var_current <- var_current[[selected]]

# getting geographic polygons of areas of interest
## simple polygons of world countries (only for plotting purposes)
wld <- world(resolution = 3, level = 0, path = "Data")

## USA Lower 48 states
us <- gadm(country = "USA", level = 1, path = "Data")
us <- us[!us$NAME_1 %in% c("Alaska", "Hawaii"), ]
us <- aggregate(us)

## Mexico
mx <- gadm(country = "MEX", level = 0, path = "Data")

## extent for a large portion of the Americas
ame_ext <- ext(-135, -35, -60, 80)
# ------------------------------------------------------------------------------


# Preparing data for analyses --------------------------------------------------
# preparing variables in areas and scenarios for comparisons
## USA vs Americas
usa <- crop(var_current, us, mask = TRUE)
americas <- crop(var_current, ame_ext)

## Ecuador current vs future
au_var <- crop(var_current, au_area, mask = TRUE)
mex <- crop(var_current, mx, mask = TRUE)

## values of variables in areas of interest
### USA vs Americas
usav <- as.data.frame(usa)
americasv <- as.data.frame(americas)

### accessible to A. ultramarina vs Mexico
au_varv <- as.data.frame(au_var)
mexv <- as.data.frame(mex)

# preparing a sample background to emulated random sampling as in Maxent (USA)
## number of pixels with data in USA
npusa <- nrow(usav)

## n = 5000, 10000, 20000
set.seed(1)
usav5 <- usav[sample(npusa, 5000, replace = FALSE), ]
usav10 <- usav[sample(npusa, 10000, replace = FALSE), ]
usav20 <- usav[sample(npusa, 20000, replace = FALSE), ]

## exploring density of points in examples
var_den <- c(2, 5)

### USA vs Americas
usak <- kde(usav[, var_den])
americask <- kde(americasv[, var_den])

### accessible to A. ultramarina vs Mexico
au_vark <- kde(au_varv[, var_den])
mexk <- kde(mexv[, var_den])

## creating blocks of points in environmental space
## preparing sets of points (small blocks in environmental space) to detect 
## differences in the way dissimilarities can be measured using options for
## variable processing and type of distance

### number of columns and rows for block matrix
n_cols <- 20
n_rows <- 20

### subset of variables used for better visualization
data <- americasv[, var_den]
id <- paste(data[, 1], data[, 2])

### determining critical values for block matrix
xrange <- range(data[, 1])
xinter <- diff(xrange)/n_cols
yrange <- range(data[, 2])
yinter <- diff(yrange)/n_rows
xlb <- seq(xrange[1], xrange[2], xinter)
xlb[length(xlb)] <- xrange[2]
ylb <- seq(yrange[1], yrange[2], yinter)
ylb[length(ylb)] <- yrange[2]

### assigning block ids according to previous specifications
blocks <- assign_blocks(data, 1, 2, n_cols, n_rows, xlb, ylb, 
                        block_type = "equal_area")
blocks <- blocks[match(id, paste(blocks[, 1], blocks[, 2])), ]

### plotting points in calibration and projection areas
ublocks <- unique(blocks[, 3])
colsb <- sample(purplow(length(ublocks)))

plot(blocks[, 1:2], col = colsb[as.factor(blocks[, 3])])
points(usav[, var_den], pch = 16, col = alpha("gray75", 0.3))

### labels of blocks produced previously
cents <- lapply(ublocks, function(x) {
  sel <- blocks[blocks[, 3] == x, 1:2]
  if (class(sel)[1] %in% c("matrix", "data.frame")) {
    cen <- apply(sel, 2, mean)
  } else {
    cen <- sel
  }
  
  text(cen[1], cen[2], labels = x, cex = 0.5)
})

### selecting blocks of interest
sel_block <- c(211, 317, 149, 215, 127, 429)

block_sel <- blocks[, 3] %in% sel_block

### corners of blocks
block_coor <- list(c(-14.6, 0.0, -11.03, 24.2), c(3.26, 24.2, 6.83, 48.4),
                   c(-25.32, 24.2, -21.75, 48.4), c(-14.6, 96.8, -11.03, 121.0),
                   c(-28.89, 0.0, -25.32, 24.2), c(21.12, 193.6, 24.70, 217.8))

### blocks of interest (centroid coordinates)
regions <- t(sapply(block_coor, function(x) {
  c(mean(x[c(1, 3)]), mean(x[c(2, 4)]))
}))
# ------------------------------------------------------------------------------


# Running analyses -------------------------------------------------------------
## We will ask the function to produce detailed information about:
##   - how many variables are found to present non-analogous conditions
##   - how non-analogous conditions are found, per variable and for multiple 
##     variables combined


# Experiment 1
## MOP analysis with six variables to demonstrate the way variables outside 
## ranges of reference conditions can be identified with more detail,
## no distances calculated
mop_ame <- mop(m = usa, g = americas, type = "detailed")
mop_mex <- mop(m = au_var, g = mex, type = "detailed")


# Experiment 2
## examples of MOP analyses using options related to type of distance used 
## and variable scaling which can generate variability in results. We used only 
## two variables to be able to show results projected into environmental space 

## Non-rescaled distance
### distance measured to the 10% closest points in the calibration area
mop_ames_eucl10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                       type = "detailed", calculate_distance = TRUE, 
                       where_distance = "all", percentage = 10, parallel = TRUE, n_cores = 4)
mop_ames_eu_s10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                       type = "detailed", calculate_distance = TRUE, 
                       where_distance = "all", percentage = 10, 
                       scale = TRUE, center = TRUE, parallel = TRUE, n_cores = 4)
mop_ames_maha10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                       type = "detailed", calculate_distance = TRUE, 
                       where_distance = "all", distance = "mahalanobis", 
                       percentage = 10, parallel = TRUE, n_cores = 4)

### distance measured to the 5% closest points in the calibration area
mop_ames_eucl5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 5, parallel = TRUE, n_cores = 4)
mop_ames_eu_s5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 5, 
                      scale = TRUE, center = TRUE, parallel = TRUE, n_cores = 4)
mop_ames_maha5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", distance = "mahalanobis",
                      percentage = 5, parallel = TRUE, n_cores = 4)

### distance measured to the 1% closest points in the calibration area
mop_ames_eucl1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 1, parallel = TRUE, n_cores = 4)
mop_ames_eu_s1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 1, 
                      scale = TRUE, center = TRUE, parallel = TRUE, n_cores = 4)
mop_ames_maha1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", distance = "mahalanobis",
                      percentage = 1, parallel = TRUE, n_cores = 4)

## MOP dissimilarity values rescaled 0-1
### distance measured to the 10% closest points in the calibration area
mop_ame_eucl10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 10, 
                      rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_eu_s10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed", calculate_distance = TRUE, 
                      where_distance = "all", percentage = 10, 
                      scale = TRUE, center = TRUE, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_maha10 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                      type = "detailed",  calculate_distance = TRUE, 
                      where_distance = "all", distance = "mahalanobis", 
                      percentage = 10, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)

### distance measured to the 5% closest points in the calibration area
mop_ame_eucl5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", percentage = 5,
                     rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_eu_s5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", percentage = 5,  
                     scale = TRUE, center = TRUE, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_maha5 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", distance = "mahalanobis", 
                     percentage = 5, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)

### distance measured to the 1% closest points in the calibration area
mop_ame_eucl1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", percentage = 1, 
                     rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_eu_s1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", percentage = 1, 
                     scale = TRUE, center = TRUE, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)
mop_ame_maha1 <- mop(m = usa[[var_den]], g = americas[[var_den]], 
                     type = "detailed", calculate_distance = TRUE, 
                     where_distance = "all", distance = "mahalanobis",
                     percentage = 1, rescale_distance = TRUE, parallel = TRUE, n_cores = 4)


# Experiment 3
## MOP analysis with six raw variables, distances to the 10% closest m points
## using sampled backgrounds and all data, USA only
mop_ame5 <- mop(m = usav5, g = americas, type = "detailed",
                calculate_distance = TRUE, percentage = 10,
                parallel = TRUE, n_cores = 4)
mop_ame10 <- mop(m = usav10, g = americas, type = "detailed", 
                 calculate_distance = TRUE, percentage = 10,
                 parallel = TRUE, n_cores = 4)
mop_ame20 <- mop(m = usav20, g = americas, type = "detailed", 
                 calculate_distance = TRUE, percentage = 10,
                 parallel = TRUE, n_cores = 4)
mop_amea <- mop(m = usa, g = americas, type = "detailed", 
                calculate_distance = TRUE, percentage = 10,
                parallel = TRUE, n_cores = 4)

# saving results
dir.create("Results")

save(mop_ame, mop_mex, mop_ame5, mop_ame10, mop_ame20, mop_amea, mop_ame_eu_s1, 
     mop_ame_eu_s5, mop_ame_eu_s10, mop_ame_eucl1, mop_ame_eucl5, 
     mop_ame_eucl10, mop_ame_maha1, mop_ame_maha5, mop_ame_maha10, 
     mop_ames_eu_s5, mop_ames_eu_s10, mop_ames_eucl1, mop_ames_eucl5, 
     mop_ames_eucl10, mop_ames_maha1, mop_ames_maha5, mop_ames_maha10, 
     file = "Results/MOPs.RData")
# ------------------------------------------------------------------------------


# Exploring results ------------------------------------------------------------
# information included in results
## a summary of the data used
mop_ame$summary

## a layer with a traditional MOP output
mop_ame$mop_basic

## a layer showing how many variables present non-analogous conditions
mop_ame$mop_simple

## several layers with detailed information of how non-analogous conditions are
## found
### non-analogous areas per variable, low end
mop_ame$mop_detailed$towards_low_end

### non-analogous areas per variable, high end
mop_ame$mop_detailed$towards_high_end

### non-analogous areas combining multiple variables, low end
mop_ame$mop_detailed$towards_low_combined

### non-analogous areas combining multiple variables, high end
mop_ame$mop_detailed$towards_high_combined
# ------------------------------------------------------------------------------


# Differences in MOP in regions of interest ------------------------------------
# getting MOP values
## non-rescaled MOPs
b_mop_values <- lapply(sel_block, function(x) {
  cbind(E_RV_10 = as.data.frame(mop_ames_eucl10$mop_distance)[blocks[, 3] == x, 1],
        E_RV_5 = as.data.frame(mop_ames_eucl5$mop_distance)[blocks[, 3] == x, 1],
        E_RV_1 = as.data.frame(mop_ames_eucl1$mop_distance)[blocks[, 3] == x, 1],
        E_SV_10 = as.data.frame(mop_ames_eu_s10$mop_distance)[blocks[, 3] == x, 1],
        E_SV_5 = as.data.frame(mop_ames_eu_s5$mop_distance)[blocks[, 3] == x, 1],
        E_SV_1 = as.data.frame(mop_ames_eu_s1$mop_distance)[blocks[, 3] == x, 1],
        M_RV_10 = as.data.frame(mop_ames_maha10$mop_distance)[blocks[, 3] == x, 1],
        M_RV_5 = as.data.frame(mop_ames_maha5$mop_distance)[blocks[, 3] == x, 1],
        M_RV_1 = as.data.frame(mop_ames_maha1$mop_distance)[blocks[, 3] == x, 1])
})

## re-scaled MOPs
b_mop_values_r <- lapply(sel_block, function(x) {
  cbind(E_RV_10 = as.data.frame(mop_ame_eucl10$mop_distance)[blocks[, 3] == x, 1],
        E_RV_5 = as.data.frame(mop_ame_eucl5$mop_distance)[blocks[, 3] == x, 1],
        E_RV_1 = as.data.frame(mop_ame_eucl1$mop_distance)[blocks[, 3] == x, 1],
        E_SV_10 = as.data.frame(mop_ame_eu_s10$mop_distance)[blocks[, 3] == x, 1],
        E_SV_5 = as.data.frame(mop_ame_eu_s5$mop_distance)[blocks[, 3] == x, 1],
        E_SV_1 = as.data.frame(mop_ame_eu_s1$mop_distance)[blocks[, 3] == x, 1],
        M_RV_10 = as.data.frame(mop_ame_maha10$mop_distance)[blocks[, 3] == x, 1],
        M_RV_5 = as.data.frame(mop_ame_maha5$mop_distance)[blocks[, 3] == x, 1],
        M_RV_1 = as.data.frame(mop_ame_maha1$mop_distance)[blocks[, 3] == x, 1])
})


# table with descriptive statistics
## non-rescaled results
stats_blocks <- sapply(b_mop_values, function(x) {
  means <- round(colMeans(x), 3)
  cis <- round(apply(x, 2, quantile, prob = c(0.025, 0.975)), 3)
  paste0(means, " (", apply(cis, 2, paste0, collapse = "-"), ")")
})

colnames(stats_blocks) <- paste0("Block ", 1:6, "; mean (CI)")
stats_blocks <- cbind(`MOP case` = colnames(b_mop_values[[1]]), stats_blocks)

## re-scaled results
stats_blocks_r <- sapply(b_mop_values_r, function(x) {
  means <- round(colMeans(x), 5)
  cis <- round(apply(x, 2, quantile, prob = c(0.025, 0.975)), 5)
  paste0(means, " (", apply(cis, 2, paste0, collapse = "-"), ")")
})

colnames(stats_blocks_r) <- paste0("Block ", 1:6, "; mean (CI)")
stats_blocks_r <- cbind(`MOP case` = colnames(b_mop_values[[1]]), stats_blocks_r)

# saving results
save(b_mop_values, b_mop_values_r, stats_blocks, stats_blocks_r, 
     file = "Results/MOP_values.RData")

# writing csv files
write.csv(stats_blocks, "Results/MOP_block_stats_nonrescaled.csv", 
          row.names = FALSE)
write.csv(stats_blocks_r, "Results/MOP_block_stats_rescaled.csv", 
          row.names = FALSE)
# ------------------------------------------------------------------------------


# Figure 2 ---------------------------------------------------------------------
# line thickness
lwdn <- 0.3
lwdh <- 0.6

# parts of figure two to be put together in Inkscape
png("Figures/Figure2a.png", width = 50, height = 50, units = "mm", res = 300)
par(mar = rep(0, 4)); plot(au_area, axes = FALSE, mar = NA, lwd = lwdh)
dev.off()

png("Figures/Figure2b.png", width = 100, height = 50, units = "mm", res = 300)
par(mar = rep(0, 4)); plot(us, axes = FALSE, mar = NA, lwd = lwdh)
dev.off()

png("Figures/Figure2c.png", width = 50, height = 50, units = "mm", res = 300)
par(mar = rep(0, 4))
plot(mx, border = NA, axes = FALSE, mar = NA, ylim = c(10, 40), lwd = lwdn) 
plot(wld, add = TRUE, border = alpha("gray65", 0.6), lwd = lwdn)
plot(mx, add = TRUE, lwd = lwdh)
plot(au_area, add = TRUE, border = "gray35", lwd = lwdh)
dev.off()

png("Figures/Figure2d.png", width = 50, height = 70, units = "mm", res = 300)
par(mar = rep(0, 4))
plot(ame_ext, mar = NA, axes = FALSE, lwd = lwdh) 
plot(wld, add = TRUE, border = alpha("gray45", 0.9), lwd = lwdn)
plot(us, add = TRUE, border = "gray15", lwd = lwdh)
dev.off()

png("Figures/Figure2e.png", width = 50, height = 50, units = "mm", res = 300)
par(mar = rep(0, 4))
plot(mop_mex$mop_simple, col = gray.colors(6)[2:5], mar = NA, axes = FALSE) 
plot(wld, add = TRUE, border = alpha("gray65", 0.6), lwd = lwdn)
plot(mx, add = TRUE, lwd = lwdh)
dev.off()

png("Figures/Figure2f.png", width = 50, height = 70, units = "mm", res = 300)
par(mar = rep(0, 4))
plot(mop_ame$mop_simple, col = gray.colors(6)[2:5], mar = NA, axes = FALSE) 
plot(wld, add = TRUE, border = alpha("gray65", 0.6), lwd = lwdn)
plot(us, add = TRUE, lwd = lwdh)
dev.off()

png("Figures/Figure2g.png", width = 50, height = 70, units = "mm", res = 300)
par(mar = rep(0, 4))
plot(mop_amea$mop_distances, col = gray.colors(250)[240:40], mar = NA, 
     axes = FALSE, legend = FALSE)
plot(mop_amea$mop_basic, col = "gray20", add = TRUE, axes = FALSE, 
     legend = FALSE) 
plot(wld, add = TRUE, border = alpha("gray65", 0.6), lwd = lwdn)
plot(us, add = TRUE, lwd = lwdh)
dev.off()
# ------------------------------------------------------------------------------


# Figure 3 ---------------------------------------------------------------------
# colors
colcal <- "#46A7F5"
colpro <- "#848687"
colsel <- "#F90616"
colcalden <- grey.colors(1000)
colcalden[1] <- NA
colproden <- grey.colors(1000)
colproden[1] <- NA
colcalden1 <- grey.colors(255)
colcalden1[1] <- NA
colproden1 <- grey.colors(255)
colproden1[1] <- NA

# exploring examples in environmental space 
png("Figures/Figure3.png", res = 600, width = 120, height = 170, units = "mm")
par(mfrow = c(3, 2), mar = c(4.2, 4.2, 0.6, 0.5), cex = 0.5)

## USA vs Americas example
limsame <- apply(americasv[, var_den], 2, range)

plot(americasv[, var_den], pch = 16, col = alpha(colcal, 0.6), 
     xlab = "", ylab = "Precipitation of driest month")
points(usav[, var_den], pch = 16, col = alpha(colpro, 0.5))
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
text(regions[, 1], regions[, 2], labels = 1:nrow(regions), cex = 0.8)
legend("topleft", cex = 1, bty = "n", 
       legend = c("Calibration area (USA)", "Transfer area (Americas)"),
       fill = c(colpro, colcal))
legend(x = -50, y = 430, cex = 0.8, bty = "n", 
       legend = c("1 Inside cloud & ranges (high density)", 
                  "2 Inside cloud & ranges (low density)",
                  "3 Outside cloud, inside ranges",
                  "   (near high density)",
                  "4 Outside cloud, inside ranges",
                  "   (near low density)",
                  "5 Outside cloud & one range",
                  "6 Outside cloud & two ranges"))

## A. u. accessible area vs MEX
lims <- apply(rbind(au_varv[, var_den], mexv[, var_den]), 2, range)

plot(mexv[, var_den], pch = 16, col = alpha(colcal, 0.7), 
     xlab = "", ylab = "",
     xlim = lims[, 1], ylim = lims[, 2])
points(au_varv[, var_den], pch = 16, col = alpha(colpro, 0.5))
legend("topleft", legend = c("Accessible area", 
                             "Transfer area (Mexico)"),
       fill = c(colpro, colcal), cex = 1, bty = "n")

## USA vs Americas (density of points)
plot(limsame, type = "n", xlab = "", ylab = "Precipitation of driest month")
image(usak$eval.points[[1]], usak$eval.points[[2]], z = usak$estimate, 
        col = colcalden, add = TRUE)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend_bar(position = "topleft", col = colcalden, title = "Density", cex = 1.1)
legend("topright", cex = 1, bty = "n", legend = "USA  ")

## Ecuador current vs future (density of points)
plot(lims, type = "n", xlab = "", ylab = "")
image(au_vark$eval.points[[1]], au_vark$eval.points[[2]], z = au_vark$estimate, 
      col = colcalden1, add = TRUE)
legend("topleft", cex = 1, bty = "n", legend = "Accessible area")

## USA vs Americas (density of points)
plot(limsame, type = "n", xlab = "Minimum temperature", 
     ylab = "Precipitation of driest month")
image(americask$eval.points[[1]], americask$eval.points[[2]], 
      z = americask$estimate, col = colproden, add = TRUE)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend_bar(position = "topleft", col = colproden, title = "Density", 
           cex = 1.1)
legend("topright", cex = 1, bty = "n", legend = "Americas  ")

## Ecuador current vs future (density of points)
plot(lims, type = "n", xlab = "Minimum temperature", ylab = "")
image(mexk$eval.points[[1]], mexk$eval.points[[2]], z = mexk$estimate, 
      col = colproden1, add = TRUE)
legend("topleft", cex = 1, bty = "n", legend = "Mexico")

dev.off()
# ------------------------------------------------------------------------------


# Figure 4 ---------------------------------------------------------------------
# margins
martext <- rep(0, 4)
marplot <- c(0.7, 4.5, 0.1, 0.1)
marplot1 <- c(4.2, 4.5, 0.1, 0.1)
marplot2 <- c(0.7, 0.7, 0.1, 0.1)
marplot3 <- c(4.2, 0.7, 0.1, 0.1)

# colors
ceuc <- raster_shared_colors(c(mop_ame_eucl10$mop_distance, 
                               mop_ame_eucl5$mop_distance, 
                               mop_ame_eucl1$mop_distance,
                               mop_ame_eu_s10$mop_distance, 
                               mop_ame_eu_s5$mop_distance, 
                               mop_ame_eu_s1$mop_distance,
                               mop_ame_maha10$mop_distance, 
                               mop_ame_maha5$mop_distance, 
                               mop_ame_maha1$mop_distance), 
                             color_palette = daright, raster_plot = FALSE)

colsel <- "#F90616"

# size proportions
cext <- 0.95
cexe <- 0.8
cexl <- 0.7
cexb <- 1.2
ptcex <- 0.5

# figure 4
png("Figures/Figure4.png", res = 600, width = 166, height = 140, units = "mm")

mlay <- matrix(1:16, nrow = 4)
layout(mlay, widths = c(1, 10, 8, 8), heights = c(1, 8, 8, 10))
par(cex = 0.55)

par(mar = martext)
plot.new()
plot.new(); text(0.5, 0.55, "10% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.55, "5% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.6, "1% reference", cex = 1.2, srt = 90)

plot.new(); text(0.6, 0.5, "Euclidean (raw variables)", cex = 1.2)

par(mar = marplot)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[1]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = TRUE)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[2]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = TRUE)

par(mar = marplot1)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[3]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = TRUE)
legend_bar(position = "topleft", col = rev(daright(255)), heigh_prop = 0.3,
           width_prop = 0.05, title = "Dissimilarity", cex = 1.2)


par(mar = martext)
plot.new(); text(0.5, 0.5, "Euclidean (scaled variables)", cex = 1.2)

par(mar = marplot2)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[4]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[5]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)

par(mar = marplot3)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[6]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = FALSE)

par(mar = martext)
plot.new(); text(0.5, 0.5, "Mahalanobis (raw variables)", cex = 1.2)

par(mar = marplot2)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[7]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[8]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)

par(mar = marplot3)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[9]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = FALSE)

dev.off()
# ------------------------------------------------------------------------------


# Figure 5 ---------------------------------------------------------------------
# margins and colors
martext <- rep(0, 4)
marplot <- c(0.7, 4.5, 0.1, 0.1)
marplot1 <- c(4.2, 4.5, 0.1, 0.1)
marplot2 <- c(0.7, 0.7, 0.1, 0.1)
marplot3 <- c(4.2, 0.7, 0.1, 0.1)

ylims <- apply(do.call(rbind, b_mop_values), 2, range)

colbox <- "gray75"

# figure 5
png("Figures/Figure5.png", res = 600, width = 166, height = 140, units = "mm")

mlay <- matrix(1:16, nrow = 4)
layout(mlay, widths = c(1, 10, 8, 8), heights = c(1, 8, 8, 10))
par(cex = 0.55)

par(mar = martext)
plot.new()
plot.new(); text(0.5, 0.55, "10% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.55, "5% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.6, "1% reference", cex = 1.2, srt = 90)

plot.new(); text(0.6, 0.5, "Euclidean (raw variables)", cex = 1.2)

par(mar = marplot)
block_bxplot(b_mop_values, column = 1, axis1 = FALSE, axis2 = TRUE, 
             ylim = ylims[, 1])
block_bxplot(b_mop_values, column = 2, axis1 = FALSE, axis2 = TRUE, 
             ylim = ylims[, 1])

par(mar = marplot1)
block_bxplot(b_mop_values, column = 3, axis1 = TRUE, axis2 = TRUE, 
             ylim = ylims[, 1])

par(mar = martext)
plot.new(); text(0.5, 0.5, "Euclidean (scaled variables)", cex = 1.2)

par(mar = marplot2)
block_bxplot(b_mop_values, column = 4, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 4])
block_bxplot(b_mop_values, column = 5, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 4])

par(mar = marplot3)
block_bxplot(b_mop_values, column = 6, axis1 = TRUE, axis2 = FALSE, 
             ylim = ylims[, 4])

par(mar = martext)
plot.new(); text(0.5, 0.5, "Mahalanobis (raw variables)", cex = 1.2)

par(mar = marplot2)
block_bxplot(b_mop_values, column = 7, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 7])
block_bxplot(b_mop_values, column = 8, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 7])

par(mar = marplot3)
block_bxplot(b_mop_values, column = 9, axis1 = TRUE, axis2 = FALSE, 
             ylim = ylims[, 7])

dev.off()
# ------------------------------------------------------------------------------


# Figure 6 ---------------------------------------------------------------------
# colors and legends
vame <- as.numeric(unlist(unique(mop_ame$mop_simple)))
vmex <- as.numeric(unlist(unique(mop_mex$mop_simple)))

ulowame <- unique(mop_ame$mop_detailed$towards_low_combined)
llowame <- gsub("_", "", ulowame[, 1])

uhighame <- unique(mop_ame$mop_detailed$towards_high_combined)
lhighame <- gsub("_", "", uhighame[, 1])

ulowau <- unique(mop_mex$mop_detailed$towards_low_combined)
llowau <- gsub("_", "", ulowau[, 1])

uhighau <- unique(mop_mex$mop_detailed$towards_high_combined)
lhighau <- gsub("_", "", uhighau[, 1])

cosim <- rev(purplow(4))
clowame <- sample(bluered(nrow(ulowame)))
chighame <- sample(bluered(nrow(uhighame)))
clowau <- sample(bluered(nrow(ulowau)))
chighau <- sample(bluered(nrow(uhighau)))

# margins
marplot <- c(0.1, 0.1, 0.1, 0.1)
martext <- rep(0, 4)

# the figure
png("Figures/Figure6.png", res = 600, width = 120, height = 185, units = "mm")
mat <- matrix(1:12, nrow = 4)
layout(mat, widths = c(1, rep(10, 3)), heights = c(1, rep(10, 3)))

## labels Y
par(mar = martext)
par(cex = 0.6)
plot.new()
plot.new() 
text(0.5, 0.5, "Number of variables (non-analogous)", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.5, "Non-analogous towards low values", cex = 1.2, 
                 srt = 90)
plot.new(); text(0.5, 0.5, "Non-analogous towards high values", cex = 1.2, 
                 srt = 90)

# USA vs Americas 
plot.new(); text(0.5, 0.5, "USA compared to Americas", cex = 1.2)

par(mar = marplot, cex = 0.8)

plot(mop_ame$mop_simple, col = cosim, axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(us, add = TRUE)
box()
legend("bottomleft", legend = vame, 
       fill = cosim, cex = 0.7, bty = "n")

plot(mop_ame$mop_detailed$towards_low_combined, col = clowame, 
     axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(us, add = TRUE)
box()
legend("bottomleft", legend = llowame, 
       fill = clowame, cex = 0.7, bty = "n")

plot(mop_ame$mop_detailed$towards_high_combined, col = chighame, 
     axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(us, add = TRUE)
box()
legend("bottomleft", legend = lhighame, 
       fill = chighame, cex = 0.7, bty = "n")

# Accessible area vs MEX 
par(mar = martext)
par(cex = 0.6)
plot.new(); text(0.5, 0.5, "Accessible area compared to MEX", cex = 1.2)

par(mar = marplot, cex = 0.8)

plot(mop_mex$mop_simple, col = cosim, axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(au_area, add = TRUE)
box()
legend("bottomleft", legend = "Reference area", 
       fill = NA, cex = 0.7, bty = "n")

plot(mop_mex$mop_detailed$towards_low_combined, col = clowau, 
     axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(au_area, add = TRUE)
box()
legend("bottomleft", legend = llowau, 
       fill = clowau, cex = 0.7, bty = "n")

plot(mop_mex$mop_detailed$towards_high_combined, col = chighau, 
     axes = F, legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
plot(au_area, add = TRUE)
box()
legend("bottomleft", legend = lhighau, 
       fill = chighau, cex = 0.7, bty = "n")

dev.off()
# ------------------------------------------------------------------------------


# Figure 7 ---------------------------------------------------------------------
# figure of classic results of MOP analyses
## margins 
marplot <- c(0.1, 0.1, 0.1, 0.1)
martext <- rep(0, 4)

## preparing colors
cbas <- raster_shared_colors(c(mop_amea$mop_distance, mop_ame20$mop_distance, 
                               mop_ame10$mop_distance, mop_ame5$mop_distance), 
                             color_palette = daright, raster_plot = TRUE)

csim <- raster_shared_colors(c(mop_amea$mop_simple, mop_ame20$mop_simple, 
                               mop_ame10$mop_simple, mop_ame5$mop_simple), 
                             color_palette = purplow, raster_plot = TRUE)

simp_leg <- unlist(terra::unique(c(mop_amea$mop_simple, mop_ame20$mop_simple, 
                                   mop_ame10$mop_simple, mop_ame5$mop_simple)))
simp_leg <- unique(na.omit(simp_leg))


## the figure
png("Figures/Figure7.png", res = 600, width = 160, height = 120, units = "mm")

mat <- matrix(1:15, nrow = 3)
layout(mat, widths = c(1, rep(10, 4)), heights = c(1, 10, 10))

## labels Y
par(mar = martext)
par(cex = 0.6)
plot.new()
plot.new(); text(0.5, 0.5, "Dissimilarity and non-analogous", cex = 1.2, 
                 srt = 90)
plot.new(); text(0.5, 0.5, "Number of variables (non-analogous)", cex = 1.2, 
                 srt = 90)

## All data
plot.new(); text(0.5, 0.5, "All data (40,207)", cex = 1.2)

par(mar = marplot)
plot(mop_amea$mop_distances, col = cbas[[1]], axes = F, legend = F, mar = NA)
image(mop_amea$mop_basic, col = "#000000", add = TRUE)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = "Non-analogous", fill = "#000000", bty = "n", 
       cex = 0.8)
legend_bar(position = c(-130, -50), col = rev(daright(255)), 
           title = "Dissimilarity", cex = 1.2)

plot(mop_amea$mop_simple, col = csim[[1]], axes = F, legend = F, mar = NA)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = simp_leg,
       fill = rev(purplow(length(simp_leg))), bty = "n", 
       cex = 0.8)

## 20000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (20,000)", cex = 1.2)
par(mar = marplot)

plot(mop_ame20$mop_distances, col = cbas[[2]], axes = F, legend = F, mar = NA)
image(mop_ame20$mop_basic, col = "#000000", add = TRUE)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

plot(mop_ame20$mop_simple, col = csim[[2]], axes = F, legend = F, mar = NA)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

## 10000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (10,000)", cex = 1.2)
par(mar = marplot)

plot(mop_ame10$mop_distances, col = cbas[[3]], axes = F, legend = F, mar = NA)
image(mop_ame10$mop_basic, col = "#000000", add = TRUE)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

plot(mop_ame10$mop_simple, col = csim[[3]], axes = F, legend = F, mar = NA)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

## 5000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (5,000)", cex = 1.2)
par(mar = marplot)

plot(mop_ame5$mop_distances, col = cbas[[4]], axes = F, legend = F, mar = NA)
image(mop_ame5$mop_basic, col = "#000000", add = TRUE)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

plot(mop_ame5$mop_simple, col = csim[[4]], axes = F, legend = F, mar = NA)
plot(wld, add = TRUE, border = alpha("gray55", 0.6))
box()

dev.off()
# ------------------------------------------------------------------------------


# Figure S1 ---------------------------------------------------------------------
# margins
martext <- rep(0, 4)
marplot <- c(0.7, 4.5, 0.1, 0.1)
marplot1 <- c(4.2, 4.5, 0.1, 0.1)
marplot2 <- c(0.7, 0.7, 0.1, 0.1)
marplot3 <- c(4.2, 0.7, 0.1, 0.1)

# colors
ceuc <- raster_shared_colors(c(mop_ames_eucl10$mop_distance, 
                               mop_ames_eucl5$mop_distance, 
                               mop_ames_eucl1$mop_distance), 
                             color_palette = daright, raster_plot = FALSE)

ceus <- raster_shared_colors(c(mop_ames_eu_s10$mop_distance, 
                               mop_ames_eu_s5$mop_distance, 
                               mop_ames_eu_s1$mop_distance), 
                             color_palette = daright, raster_plot = FALSE)

cmah <- raster_shared_colors(c(mop_ames_maha10$mop_distance, 
                               mop_ames_maha5$mop_distance, 
                               mop_ames_maha1$mop_distance), 
                             color_palette = daright, raster_plot = FALSE)

# size proportions
cext <- 0.95
cexe <- 0.8
cexl <- 0.7
cexb <- 1.2
ptcex <- 0.5

# figure 5
png("Figures/FigureS1.png", res = 600, width = 166, height = 140, units = "mm")

mlay <- matrix(1:16, nrow = 4)
layout(mlay, widths = c(1, 10, 8, 8), heights = c(1, 8, 8, 10))
par(cex = 0.55)

par(mar = martext)
plot.new()
plot.new(); text(0.5, 0.55, "10% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.55, "5% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.6, "1% reference", cex = 1.2, srt = 90)

plot.new(); text(0.6, 0.5, "Euclidean (raw variables)", cex = 1.2)

par(mar = marplot)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[1]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = TRUE)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[2]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = TRUE)

par(mar = marplot1)
block_scatter(americasv[, var_den], block_coor, ptcol = ceuc[[3]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = TRUE)
legend_bar(position = "topleft", col = rev(daright(255)), heigh_prop = 0.3,
           width_prop = 0.05, title = "Dissimilarity", cex = 1.2)

par(mar = martext)
plot.new(); text(0.5, 0.5, "Euclidean (scaled variables)", cex = 1.2)

par(mar = marplot2)
block_scatter(americasv[, var_den], block_coor, ptcol = ceus[[1]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)
block_scatter(americasv[, var_den], block_coor, ptcol = ceus[[2]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)

par(mar = marplot3)
block_scatter(americasv[, var_den], block_coor, ptcol = ceus[[3]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = FALSE)

par(mar = martext)
plot.new(); text(0.5, 0.5, "Mahalanobis (raw variables)", cex = 1.2)

par(mar = marplot2)
block_scatter(americasv[, var_den], block_coor, ptcol = cmah[[1]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)
block_scatter(americasv[, var_den], block_coor, ptcol = cmah[[2]], 
              ptcex = ptcex, colsel = colsel, axis1 = FALSE, axis2 = FALSE)

par(mar = marplot3)
block_scatter(americasv[, var_den], block_coor, ptcol = cmah[[3]], 
              ptcex = ptcex, colsel = colsel, axis1 = TRUE, axis2 = FALSE)

dev.off()
# ------------------------------------------------------------------------------


# Figure S2 ---------------------------------------------------------------------
# margins and colors
martext <- rep(0, 4)
marplot <- c(0.7, 4.5, 0.1, 0.1)
marplot1 <- c(4.2, 4.5, 0.1, 0.1)
marplot2 <- c(0.7, 0.7, 0.1, 0.1)
marplot3 <- c(4.2, 0.7, 0.1, 0.1)

ylims <- apply(do.call(rbind, b_mop_values_r), 2, range)

colbox <- "gray75"

# figure S2
png("Figures/FigureS2.png", res = 600, width = 166, height = 140, units = "mm")

mlay <- matrix(1:16, nrow = 4)
layout(mlay, widths = c(1, 10, 8, 8), heights = c(1, 8, 8, 10))
par(cex = 0.55)

par(mar = martext)
plot.new()
plot.new(); text(0.5, 0.55, "10% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.55, "5% reference", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.6, "1% reference", cex = 1.2, srt = 90)

plot.new(); text(0.6, 0.5, "Euclidean (raw variables)", cex = 1.2)

par(mar = marplot)
block_bxplot(b_mop_values_r, column = 1, axis1 = FALSE, axis2 = TRUE, 
             ylim = ylims[, 1])
block_bxplot(b_mop_values_r, column = 2, axis1 = FALSE, axis2 = TRUE, 
             ylim = ylims[, 1])

par(mar = marplot1)
block_bxplot(b_mop_values_r, column = 3, axis1 = TRUE, axis2 = TRUE, 
             ylim = ylims[, 1])

par(mar = martext)
plot.new(); text(0.5, 0.5, "Euclidean (scaled variables)", cex = 1.2)

par(mar = marplot2)
block_bxplot(b_mop_values_r, column = 4, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 4])
block_bxplot(b_mop_values_r, column = 5, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 4])

par(mar = marplot3)
block_bxplot(b_mop_values_r, column = 6, axis1 = TRUE, axis2 = FALSE, 
             ylim = ylims[, 4])

par(mar = martext)
plot.new(); text(0.5, 0.5, "Mahalanobis (raw variables)", cex = 1.2)

par(mar = marplot2)
block_bxplot(b_mop_values_r, column = 7, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 7])
block_bxplot(b_mop_values_r, column = 8, axis1 = FALSE, axis2 = FALSE, 
             ylim = ylims[, 7])

par(mar = marplot3)
block_bxplot(b_mop_values_r, column = 9, axis1 = TRUE, axis2 = FALSE, 
             ylim = ylims[, 7])

dev.off()
# ------------------------------------------------------------------------------


# Figure S3 --------------------------------------------------------------------
# margins
martext <- rep(0, 4)
marplot <- c(0.1, 0.1, 0.1, 0.1)

# variables involved
varsn <- names(americas)
varsn1 <- gsub("_", "", varsn)

# color
col0 <- "#740211"

png("Figures/FigureS3.png", res = 600, width = 166, height = 130, units = "mm")
mat <- matrix(1:35, nrow = 5, byrow = TRUE)
layout(mat, widths = c(0.9, rep(10, 6)), heights = c(0.9, rep(10, 4)))

## labes Y
par(mar = martext)
par(cex = 0.7)
plot.new()
for (i in varsn1) {
  plot.new(); text(0.5, 0.4, i, cex = 1)
}

## row1
plot.new(); text(0.6, 0.5, "Towards low values", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(mop_ame$mop_detailed$towards_low_end[[i]], col = col0, axes = F, 
       legend = F, mar = NA)
  plot(wld, border = alpha("gray55", 0.6), add = TRUE)
  box()
}

## row2
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards high values", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(mop_ame$mop_detailed$towards_high_end[[i]], col = col0, axes = F,
       legend = F, mar = NA)
  plot(wld, border = alpha("gray55", 0.6), add = TRUE)
  box()
}

## row3
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards low values", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(mop_mex$mop_detailed$towards_low_end[[i]], col = col0, axes = F, 
       legend = F, mar = NA)
  plot(wld, border = alpha("gray55", 0.6), add = TRUE)
  plot(au_area, add = TRUE)
  box()
}

## row4
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards high values", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(mop_mex$mop_detailed$towards_high_end[[i]], col = col0, 
       axes = F, legend = F, mar = NA)
  plot(wld, border = alpha("gray55", 0.6), add = TRUE)
  plot(au_area, add = TRUE)
  box()
}
legend("bottomleft", legend = "Non-analogous", fill = col0, bty = "n", 
       cex = 0.7)

dev.off()
# ------------------------------------------------------------------------------


# Figure S4 --------------------------------------------------------------------
# figure of combination of variables outside ranges
## margins 
marplot <- c(0.1, 0.1, 0.1, 0.1)
martext <- rep(0, 4)

## preparing colors
ulowa <- unique(mop_amea$mop_detailed$towards_low_combined)
llowa <- gsub("_", "", ulowa[, 1])
uhigha <- unique(mop_amea$mop_detailed$towards_high_combined)
lhigha <- gsub("_", "", uhigha[, 1])

ulow20 <- unique(mop_ame20$mop_detailed$towards_low_combined)
llow20 <- gsub("_", "", ulow20[, 1])
uhigh20 <- unique(mop_ame20$mop_detailed$towards_high_combined)
lhigh20 <- gsub("_", "", uhigh20[, 1])

ulow10 <- unique(mop_ame10$mop_detailed$towards_low_combined)
llow10 <- gsub("_", "", ulow10[, 1])
uhigh10 <- unique(mop_ame10$mop_detailed$towards_high_combined)
lhigh10 <- gsub("_", "", uhigh10[, 1])

ulow5 <- unique(mop_ame5$mop_detailed$towards_low_combined)
llow5 <- gsub("_", "", ulow5[, 1])
uhigh5 <- unique(mop_ame5$mop_detailed$towards_high_combined)
lhigh5 <- gsub("_", "", uhigh5[, 1])

clowa <- bluered(nrow(ulowa))
chigha <- bluered(nrow(uhigha))

clow20 <- bluered(nrow(ulow20))
chigh20 <- bluered(nrow(uhigh20))

clow10 <- bluered(nrow(ulow10))
chigh10 <- bluered(nrow(uhigh10))

clow5 <- bluered(nrow(ulow5))
chigh5 <- bluered(nrow(uhigh5))

## the figure
png("Figures/FigureS4.png", res = 600, width = 160, height = 120, units = "mm")

mat <- matrix(1:15, nrow = 3)
layout(mat, widths = c(1, rep(10, 4)), heights = c(1, 10, 10))

## labels Y
par(mar = martext)
par(cex = 0.6)
plot.new()
plot.new(); text(0.5, 0.5, "Non-analogous towards low values", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.5, "Non-analogous towards high values", cex = 1.2, srt = 90)

## All data
plot.new(); text(0.5, 0.5, "All data (40,207)", cex = 1.2)

par(mar = marplot)
plot(mop_amea$mop_detailed$towards_low_combined, col = clowa, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = llowa, fill = clowa, bty = "n", cex = 0.8)

plot(mop_amea$mop_detailed$towards_high_combined, col = chigha, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = lhigha, fill = chigha, bty = "n", cex = 0.8)

## 20000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (20,000)", cex = 1.2)

par(mar = marplot)
plot(mop_ame20$mop_detailed$towards_low_combined, col = clow20, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = llow20, fill = clow20, bty = "n", cex = 0.8)

plot(mop_ame20$mop_detailed$towards_high_combined, col = chigh20, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = lhigh20, fill = chigh20, bty = "n", cex = 0.8)

## 10000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (10,000)", cex = 1.2)
par(mar = marplot)

par(mar = marplot)
plot(mop_ame10$mop_detailed$towards_low_combined, col = clow10, axes = F, 
      legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = llow10, fill = clow10, bty = "n", cex = 0.8)

plot(mop_ame10$mop_detailed$towards_high_combined, col = chigh10, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = lhigh10, fill = chigh10, bty = "n", cex = 0.8)

## 5000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (5,000)", cex = 1.2)
par(mar = marplot)

par(mar = marplot)
plot(mop_ame5$mop_detailed$towards_low_combined, col = clow5, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = llow5, fill = clow5, bty = "n", cex = 0.8)

plot(mop_ame5$mop_detailed$towards_high_combined, col = chigh5, axes = F, 
     legend = F, mar = NA)
plot(wld, border = alpha("gray55", 0.6), add = TRUE)
box()
legend("bottomleft", legend = lhigh5, fill = chigh5, bty = "n", cex = 0.8)

dev.off()
# ------------------------------------------------------------------------------
