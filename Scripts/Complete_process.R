# ------------------------------------------------------------------------------
# Project: Expanding interpretation of results from the MOP metric
# Author: Marlon E. Cobos, Hannah Owens, Jorge Soberon, and A. Townsend Peterson
# Date: 19/09/2022
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
# Note: some results will be written in your working directory.
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# these lines load packages. If needed, install packages using
# install.packages("package_name")
library(mop)
library(scales)
library(ks)
#library(maps) ?
library(geodata)
# ------------------------------------------------------------------------------


# Working directory ------------------------------------------------------------
setwd("YOUR/DIRECTORY")
# ------------------------------------------------------------------------------


# Example data -----------------------------------------------------------------
# accessible area simulated for the species Aphelocoma ultramarina
au_area <- vect("Data/M_simulated.gpkg")

# downloading environmental variables 
var_current <- worldclim_global(var = "bio", res = 10, path = "Data")

## renaming variables 
names(var_current) <- gsub("wc2.1_10m_", "", names(var_current))

## selecting a subset of variables
selected <- c(5:7, 13:15)

var_current <- var_current[[selected]]

# getting geographic polygons of areas of interest
## USA Lower 48 states
us <- gadm(country = "USA", level = 1, path = "Data")
us <- us[!us$NAME_1 %in% c("Alaska", "Hawaii"), ]

## Mexico
mx <- gadm(country = "MEX", level = 0, path = "Data")

## extent for a large portion of the Americas
ame_ext <- extent(-135, -35, -60, 80)
# ------------------------------------------------------------------------------


# Preparing data for analyses --------------------------------------------------
# preparing variables in areas and scenarios for comparisons
## USA vs Americas
usa <- mask(crop(var_current, us), us)
americas <- crop(var_current, ame_ext)

## Ecuador current vs future
ecu_current <- mask(crop(var_current, ec), ec)
ecu_future <- mask(crop(var_future, ec), ec)

## values of variables in areas of interest
### USA vs Americas
usav <- na.omit(usa[])
americasv <- na.omit(americas[])

### Ecuador current vs future
ecucv <- na.omit(ecu_current[])
ecufv <- na.omit(ecu_future[])

# preparing a sample background to emulated random sampling as in Maxent (USA)
## number of pixels with data in USA
npusa <- nrow(usav)

## n = 5000, 10000, 20000
usav5 <- usav[sample(npusa, 5000, replace = FALSE), ]
usav10 <- usav[sample(npusa, 10000, replace = FALSE), ]
usav20 <- usav[sample(npusa, 20000, replace = FALSE), ]

## exploring density of points in examples
### USA vs Americas
usak <- kde(usav[, c(2, 4)])
americask <- kde(americasv[, c(2, 4)])

### Ecuador current vs future
ecuck <- kde(ecucv[, c(2, 4)])
ecufk <- kde(ecufv[, c(2, 4)])

## creating blocks of points in environmental space
## preparing sets of points (small blocks in environmental space) to detect 
## differences in the way dissimilarities can be measured using options for
## variable processing and type of distance

### number of columns and rows for block matrix
n_cols <- 20
n_rows <- 20

### subset of variables used for better visualization
data <- americasv[, c(2, 4)]
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
points(usav[, c(2, 4)], pch = 16, col = alpha("gray75", 0.3))

### labels of blocks produced previously
cents <- lapply(ublocks, function(x) {
  sel <- blocks[blocks[, 3] == x, 1:2]
  if (class(sel)[1] == "matrix") {
    cen <- apply(sel, 2, mean)
  } else {
    cen <- sel
  }
  
  text(cen[1], cen[2], labels = x, cex = 0.5)
})

### selecting blocks of interest
sel_block <- c(234, 281, 151, 343, 1, 437)

block_sel <- blocks[, 3] %in% sel_block

### corners of blocks
block_coor <- list(c(-110, 98, -74.8, 147), c(-39.6, 343, -4.4, 392),
                   c(-250.8, 147, -215.6, 196), c(66, 294, 101.2, 343),
                   c(-462, 0, -426.8, 49), c(206.8, 784, 242, 833))

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
## examples of MOP analyses using options related to type of distance used 
## and variable scaling which can generate variability in results. We used only 
## two variables to be able to show results projected into environmental space 

## Non-rescaled MOP (dissimilarity is raw distance) non-analogous areas are 10%
## more distant than other areas
### distance measured to the 5% closest points in the calibration area
mop_ames_eucl10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                       mop_type = "detailed", percent = 10, rescale_mop = FALSE)
mop_ames_eu_s10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                       mop_type = "detailed", percent = 10, rescale_mop = FALSE, 
                       scale = TRUE, center = TRUE)
mop_ames_maha10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                       mop_type = "detailed", distance = "mahalanobis", 
                       percent = 10, rescale_mop = FALSE)

### distance measured to the 5% closest points in the calibration area
mop_ames_eucl5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 5, rescale_mop = FALSE)
mop_ames_eu_s5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 5, rescale_mop = FALSE, 
                      scale = TRUE, center = TRUE)
mop_ames_maha5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", distance = "mahalanobis", 
                      percent = 5, rescale_mop = FALSE)

### distance measured to the 1% closest points in the calibration area
mop_ames_eucl1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 1, rescale_mop = FALSE)
mop_ames_eu_s1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 1, rescale_mop = FALSE, 
                      scale = TRUE, center = TRUE)
mop_ames_maha1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", distance = "mahalanobis", 
                      percent = 1, rescale_mop = FALSE)

## MOP dissimilarity values rescaled 0-1
### distance measured to the 5% closest points in the calibration area
mop_ame_eucl10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 10)
mop_ame_eu_s10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", percent = 10, 
                      scale = TRUE, center = TRUE)
mop_ame_maha10 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                      mop_type = "detailed", distance = "mahalanobis", 
                      percent = 10)

### distance measured to the 5% closest points in the calibration area
mop_ame_eucl5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", percent = 5)
mop_ame_eu_s5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", percent = 5, 
                     scale = TRUE, center = TRUE)
mop_ame_maha5 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", distance = "mahalanobis", 
                     percent = 5)

### distance measured to the 1% closest points in the calibration area
mop_ame_eucl1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", percent = 1)
mop_ame_eu_s1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", percent = 1, 
                     scale = TRUE, center = TRUE)
mop_ame_maha1 <- mop(m = usa[[c(2, 4)]], g = americas[[c(2, 4)]], 
                     mop_type = "detailed", distance = "mahalanobis", 
                     percent = 1)


# Experiment 2
## MOP analysis with six raw variables, using Euclidean distances (rescaled 0-1)
## using raster layers representing the entire USA and Ecuador
mop_ame <- mop(m = usa, g = americas, mop_type = "detailed")
mop_ecu <- mop(m = ecu_current, g = ecu_future, mop_type = "detailed")


# Experiment 3
## MOP analysis with six raw variables, using Euclidean distances (non-rescaled),
## using sampled backgrounds and all data, USA only
mop_ame5 <- mop(m = usav5, g = americas, mop_type = "detailed",
                scale = TRUE, center = TRUE, rescale_mop = FALSE)
mop_ame10 <- mop(m = usav10, g = americas, mop_type = "detailed", 
                 scale = TRUE, center = TRUE, rescale_mop = FALSE)
mop_ame20 <- mop(m = usav20, g = americas, mop_type = "detailed", 
                 scale = TRUE, center = TRUE, rescale_mop = FALSE)
mop_amea <- mop(m = usa, g = americas, mop_type = "detailed", 
                scale = TRUE, center = TRUE, rescale_mop = FALSE)

# saving results
dir.create("Results")

save(mop_ame, mop_ecu, mop_ame5, mop_ame10, mop_ame20, mop_amea, mop_ame_eu_s1, 
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
### a table to help interpret values in layers that combine multiple variables
mop_ame$mop_detailed$interpretation_combined

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
  cbind(E_RV_10 = na.omit(mop_ames_eucl10$mop_basic[])[blocks[, 3] == x],
        E_RV_5 = na.omit(mop_ames_eucl5$mop_basic[])[blocks[, 3] == x],
        E_RV_1 = na.omit(mop_ames_eucl1$mop_basic[])[blocks[, 3] == x],
        E_SV_10 = na.omit(mop_ames_eu_s10$mop_basic[])[blocks[, 3] == x],
        E_SV_5 = na.omit(mop_ames_eu_s5$mop_basic[])[blocks[, 3] == x],
        E_SV_1 = na.omit(mop_ames_eu_s1$mop_basic[])[blocks[, 3] == x],
        M_RV_10 = na.omit(mop_ames_maha10$mop_basic[])[blocks[, 3] == x],
        M_RV_5 = na.omit(mop_ames_maha5$mop_basic[])[blocks[, 3] == x],
        M_RV_1 = na.omit(mop_ames_maha1$mop_basic[])[blocks[, 3] == x])
})

## re-scaled MOPs
b_mop_values_r <- lapply(sel_block, function(x) {
  cbind(E_RV_10 = na.omit(mop_ame_eucl10$mop_basic[])[blocks[, 3] == x],
        E_RV_5 = na.omit(mop_ame_eucl5$mop_basic[])[blocks[, 3] == x],
        E_RV_1 = na.omit(mop_ame_eucl1$mop_basic[])[blocks[, 3] == x],
        E_SV_10 = na.omit(mop_ame_eu_s10$mop_basic[])[blocks[, 3] == x],
        E_SV_5 = na.omit(mop_ame_eu_s5$mop_basic[])[blocks[, 3] == x],
        E_SV_1 = na.omit(mop_ame_eu_s1$mop_basic[])[blocks[, 3] == x],
        M_RV_10 = na.omit(mop_ame_maha10$mop_basic[])[blocks[, 3] == x],
        M_RV_5 = na.omit(mop_ame_maha5$mop_basic[])[blocks[, 3] == x],
        M_RV_1 = na.omit(mop_ame_maha1$mop_basic[])[blocks[, 3] == x])
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
  means <- round(colMeans(x), 4)
  cis <- round(apply(x, 2, quantile, prob = c(0.025, 0.975)), 4)
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
cord_amex <- c(bbox(ame_ext))

map(col = NA, xlim = c(-170, -16), ylim = c(-60, 90))
image(var_current$bio7, add = TRUE, col = rev(terrain.colors(255)))
map(col = "gray65", add = TRUE)
plot(us, border = "red", lwd = 2, add = TRUE)

scalebar(d = 5000)

axis(1, at = seq(-160, -20, 20))
axis(2, las = 1)
box()

map(col = NA, xlim = c(-170, -16), ylim = c(-60, 90))
image(var_current$bio7, add = TRUE, col = rev(terrain.colors(255)))
map(col = "gray65", add = TRUE)
rect(cord_amex[1], cord_amex[2], cord_amex[3], cord_amex[4], 
     border = "blue", lwd = 2)

# ------------------------------------------------------------------------------



# Figure 3 ---------------------------------------------------------------------
# colors
colcal <- "#46A7F5"
colpro <- "#848687"
colsel <- "#F90616"
colcalden <- grey.colors(1000)
colcalden[1] <- NA
colproden <- purplow(1000)
colproden[1] <- NA
colcalden1 <- grey.colors(255)
colcalden1[1] <- NA
colproden1 <- purplow(255)
colproden1[1] <- NA

# exploring examples in environmental space 
png("Figures/Figure3.png", res = 600, width = 120, height = 170, units = "mm")
par(mfrow = c(3, 2), mar = c(4.2, 4.2, 0.6, 0.5), cex = 0.5)

## USA vs Americas example
limsame <- apply(americasv[, c(2, 4)], 2, range)

plot(americasv[, c(2, 4)], pch = 16, col = alpha(colcal, 0.6), 
     xlab = "", ylab = "Precipitation of wettest month")
points(usav[, c(2, 4)], pch = 16, col = alpha(colpro, 0.5))
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
text(regions[, 1], regions[, 2], labels = 1:nrow(regions), cex = 0.8)
legend("topleft", cex = 1, bty = "n", 
       legend = c("Calibration area (USA)", "Transfer area (Americas)"),
       fill = c(colpro, colcal))
legend(x = -510, y = 870, cex = 0.8, bty = "n", 
       legend = c("1 Inside cloud & ranges (high density)", 
                  "2 Inside cloud & ranges (low density)",
                  "3 Outside cloud, inside ranges",
                  "   (near high density)",
                  "4 Outside cloud, inside ranges",
                  "   (near low density)",
                  "5 Outside cloud & one range",
                  "6 Outside cloud & two ranges"))

## Ecuador current vs future example
lims <- apply(rbind(ecucv[, c(2, 4)], ecufv[, c(2, 4)]), 2, range)

plot(ecufv[, c(2, 4)], pch = 16, col = alpha(colcal, 0.7), 
     xlab = "", ylab = "",
     xlim = lims[, 1], ylim = lims[, 2])
points(ecucv[, c(2, 4)], pch = 16, col = alpha(colpro, 0.5))
legend("topleft", legend = c("Calibration area (ECU current)", 
                             "Transfer area (ECU future)"),
       fill = c(colpro, colcal), cex = 1, bty = "n")

## USA vs Americas (density of points)
plot(limsame, type = "n", xlab = "", ylab = "Precipitation of wettest month")
image(usak$eval.points[[1]], usak$eval.points[[2]], z = usak$estimate, 
        col = colcalden, add = TRUE)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend_bar(position = "topleft", col = colcalden, title = "Density", cex = 1.1)
legend("topright", cex = 1, bty = "n", legend = "USA  ")

## Ecuador current vs future (density of points)
plot(lims, type = "n", xlab = "", ylab = "")
image(ecuck$eval.points[[1]], ecuck$eval.points[[2]], z = ecuck$estimate, 
      col = colcalden1, add = TRUE)
legend("topleft", cex = 1, bty = "n", legend = "ECU current")

## USA vs Americas (density of points)
plot(limsame, type = "n", xlab = "Minimum temperature", 
     ylab = "Precipitation of wettest month")
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
image(ecufk$eval.points[[1]], ecufk$eval.points[[2]], z = ecufk$estimate, 
      col = colproden1, add = TRUE)
legend("topleft", cex = 1, bty = "n", legend = "ECU future")

dev.off()
# ------------------------------------------------------------------------------



# Figure 4 ---------------------------------------------------------------------
# colors
reura10 <- plot_raster_inf(mop_ames_eucl10, result = "basic", 
                           color_palette = daright)

reusc10 <- plot_raster_inf(mop_ames_eu_s10, result = "basic", 
                           color_palette = daright)

rmaha10 <- plot_raster_inf(mop_ames_maha10, result = "basic", 
                           color_palette = daright)

reura5 <- plot_raster_inf(mop_ames_eucl5, result = "basic", 
                          color_palette = daright)

reusc5 <- plot_raster_inf(mop_ames_eu_s5, result = "basic", 
                          color_palette = daright)

rmaha5 <- plot_raster_inf(mop_ames_maha5, result = "basic", 
                          color_palette = daright)

reura1 <- plot_raster_inf(mop_ames_eucl1, result = "basic", 
                          color_palette = daright)

reusc1 <- plot_raster_inf(mop_ames_eu_s1, result = "basic", 
                          color_palette = daright)

rmaha1 <- plot_raster_inf(mop_ames_maha1, result = "basic", 
                          color_palette = daright)

# size proportions
cext <- 0.95
cexe <- 0.8
cexl <- 0.7
cexb <- 1.2
ptcex <- 0.5

# initial figure showing example
png("Figures/Figure4.png", res = 600, width = 166, height = 145, units = "mm")
par(mfrow = c(3, 3), mar = c(4.2, 4.2, 0.65, 0.5), cex = 0.5)

## 10%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = reura10$col[reura10$val_factor], 
     cex = ptcex, xlab = "", ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Euclidean-raw variables)", bty = "n", 
       cex = cext)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = reusc10$col[reusc10$val_factor], 
     cex = ptcex, xlab = "", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = rmaha10$col[rmaha10$val_factor],  
     xlab = "", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

## 5%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = reura5$col[reura5$val_factor], 
     cex = ptcex, xlab = "", ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Euclidean-raw variables)", bty = "n", 
       cex = cext)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = reusc5$col[reusc5$val_factor], 
     cex = ptcex, xlab = "", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = rmaha5$col[rmaha5$val_factor],  
     xlab = "", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

## 1%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = reura1$col[reura1$val_factor], 
     cex = ptcex, xlab = "Minimum temperature", 
     ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Euclidean-raw variables)", bty = "n", 
       cex = cext)
legend(x = -480, y = 850, legend = "Non-analogous", fill = "#010101", bty = "n", 
       cex = cexe)
legend_bar(position = c(-450, 450), col = reura1$col, title = "Dissimilarity", 
           cex = cexb)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = reusc1$col[reusc1$val_factor], 
     cex = ptcex, xlab = "Minimum temperature", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = rmaha1$col[rmaha1$val_factor],  
     xlab = "Minimum temperature", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

dev.off()
# ------------------------------------------------------------------------------



# Figure 5 ---------------------------------------------------------------------
# margins
marplot <- c(4.7, 4.2, 1.5, 0.1)
martext <- rep(0, 4)
colbox <- alpha(c("#1b9e77", "#d95f02", "#7570b3"), 0.7)

# complement for figure 5
png("Figures/Figure5.png", res = 600, width = 166, height = 95, units = "mm")

mat <- matrix(1:12, nrow = 3)
layout(mat, widths = rep(10, 4), heights = c(0.7, rep(10, 2)))
par(cex = 0.45)

## boxplots of numerical results in blocks
par(mar = martext)
plot.new(); text(0.6, 0.3, "Block 1", cex = 1.2)

par(mar = marplot)
boxplot(b_mop_values[[1]], las = 2, ylab = "Distance", 
        xlab = "", outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
legend("topright", legend = paste(c("10%", "5%", "1%"), "M reference"), 
       fill = colbox, bty = "n")
title(main = "Inside cloud & ranges (high density)", font.main = 1, 
      cex.main = 0.8)

boxplot(b_mop_values_r[[1]], las = 2, ylab = "Re-scaled distance", 
        xlab = "", outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))

par(mar = martext)
plot.new(); text(0.6, 0.3, "Block 2", cex = 1.2)

par(mar = marplot)
boxplot(b_mop_values[[2]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
title(main = "Inside cloud & ranges (low density)", font.main = 1, 
      cex.main = 0.8)

boxplot(b_mop_values_r[[2]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))

par(mar = martext)
plot.new(); text(0.6, 0.3, "Block 3", cex = 1.2)

par(mar = marplot)
boxplot(b_mop_values[[3]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox, 
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
title(main = "Outside cloud, in ranges (near high density)", font.main = 1,
      cex.main = 0.8)

boxplot(b_mop_values_r[[3]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))

par(mar = martext)
plot.new(); text(0.6, 0.3, "Block 4", cex = 1.2)

par(mar = marplot)
boxplot(b_mop_values[[4]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
title(main = "Outside cloud, in ranges (near low density)", font.main = 1,
      cex.main = 0.8)

boxplot(b_mop_values_r[[4]], las = 2, ylab = "", xlab = "", 
        outline = FALSE, lwd = 0.5, col = colbox,
        pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))

dev.off()
# ------------------------------------------------------------------------------



# Figure 6 ---------------------------------------------------------------------
# colors and legends
mopame_sim <- plot_raster_inf(mop_ame, result = "simple", 
                              color_palette = purplow, reverse = TRUE)

mopecu_sim <- plot_raster_inf(mop_ecu, result = "simple", 
                              color_palette = purplow, reverse = TRUE)

mopame_tlc <- plot_raster_inf(mop_ame, result = "towards_low_combined", 
                              color_palette = bluered, reverse = FALSE)

mopame_thc <- plot_raster_inf(mop_ame, result = "towards_high_combined", 
                              color_palette = bluered, reverse = FALSE)

mopecu_tlc <- plot_raster_inf(mop_ecu, result = "towards_low_combined", 
                              color_palette = bluered, reverse = FALSE)

mopecu_thc <- plot_raster_inf(mop_ecu, result = "towards_high_combined", 
                              color_palette = bluered, reverse = FALSE)

## regions for plotting
boxpam <- t(bbox(americas@extent))
boxpamame <- SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                    proj4string = americas@crs)

boxpam <- t(bbox(ecu_current@extent))
boxpamecu <- SpatialPointsDataFrame(boxpam, data.frame(boxpam),
                                    proj4string = ecu_current@crs)

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
plot.new(); text(0.5, 0.5, "Non-analogous towards low end", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.5, "Non-analogous towards high end", cex = 1.2, srt = 90)

# USA vs Americas 
plot.new(); text(0.5, 0.5, "USA compared to Americas", cex = 1.2)

par(mar = marplot, cex = 0.8)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame$mop_simple, col = mopame_sim$col, 
      breaks = mopame_sim$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", title = "Variables", legend = mopame_sim$legend_text, 
       fill = mopame_sim$legend_col, cex = 0.7, bty = "n")

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame$mop_detailed$towards_low_combined, col = mopame_tlc$col, 
      breaks = mopame_tlc$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", title = "Variables", legend = mopame_tlc$legend_text, 
       fill = mopame_tlc$legend_col, cex = 0.7, bty = "n")

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame$mop_detailed$towards_high_combined, col = mopame_thc$col, 
      breaks = mopame_thc$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", title = "Variables", legend = mopame_thc$legend_text, 
       fill = mopame_thc$legend_col, cex = 0.7, bty = "n")

# ECU current vs future 
par(mar = martext)
par(cex = 0.6)
plot.new(); text(0.5, 0.5, "ECU current compared to future", cex = 1.2)

par(mar = marplot, cex = 0.8)

plot(boxpamecu, col = NA, axes = FALSE)
image(mop_ecu$mop_simple, col = mopecu_sim$col, 
      breaks = mopecu_sim$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomright", title = "Variables", legend = mopecu_sim$legend_text, 
       fill = mopecu_sim$legend_col, cex = 0.7, bty = "n")

plot(boxpamecu, col = NA, axes = FALSE)
image(mop_ecu$mop_detailed$towards_low_combined, col = mopecu_tlc$col, 
      breaks = mopecu_tlc$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomright", title = "Variables", legend = mopecu_tlc$legend_text, 
       fill = mopecu_tlc$legend_col, cex = 0.7, bty = "n")

plot(boxpamecu, col = NA, axes = FALSE)
image(mop_ecu$mop_detailed$towards_high_combined, col = mopecu_thc$col, 
      breaks = mopecu_thc$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomright", title = "Variables", legend = mopecu_thc$legend_text, 
       fill = mopecu_thc$legend_col, cex = 0.7, bty = "n")

dev.off()
# ------------------------------------------------------------------------------



# Figure 7 ---------------------------------------------------------------------
# figure of classic results of MOP analyses
## margins 
marplot <- c(0.1, 0.1, 0.1, 0.1)
martext <- rep(0, 4)

## preparing colors
mopame_bas <- plot_raster_inf(mop_amea, result = "basic", 
                              color_palette = daright)

mopame_bas5 <- plot_raster_inf(mop_ame5, result = "basic", 
                               color_palette = daright)

mopame_bas10 <- plot_raster_inf(mop_ame10, result = "basic", 
                                color_palette = daright)

mopame_bas20 <- plot_raster_inf(mop_ame20, result = "basic", 
                                color_palette = daright)

mopame_sim <- plot_raster_inf(mop_amea, result = "simple", 
                              color_palette = purplow, reverse = TRUE)

mopame_sim5 <- plot_raster_inf(mop_ame5, result = "simple", 
                               color_palette = purplow, reverse = TRUE)

mopame_sim10 <- plot_raster_inf(mop_ame10, result = "simple", 
                                color_palette = purplow, reverse = TRUE)

mopame_sim20 <- plot_raster_inf(mop_ame20, result = "simple", 
                                color_palette = purplow, reverse = TRUE)

## the figure
png("Figures/Figure7.png", res = 600, width = 160, height = 120, units = "mm")

mat <- matrix(1:15, nrow = 3)
layout(mat, widths = c(1, rep(10, 4)), heights = c(1, 10, 10))

## labels Y
par(mar = martext)
par(cex = 0.6)
plot.new()
plot.new(); text(0.5, 0.5, "Dissimilarity and non-analogous", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.5, "Number of variables (non-analogous)", cex = 1.2, srt = 90)

## All data
plot.new(); text(0.5, 0.5, "All data (40,207)", cex = 1.2)

par(mar = marplot)
plot(boxpamame, col = NA, axes = FALSE)
image(mop_amea$mop_basic, col = mopame_bas$col, 
      breaks = mopame_bas$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = "Non-analogous", fill = "#000000", bty = "n", 
       cex = 0.8)
legend_bar(position = c(-133, -50), col = mopame_bas$col, 
           title = "Dissimilarity", cex = 1.2)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_amea$mop_simple, col = mopame_sim$col, 
      breaks = mopame_sim$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", title = "Variables", legend = mopame_sim$legend_text,
       fill = mopame_sim$legend_col, bty = "n", 
       cex = 0.8)

## 20000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (20,000)", cex = 1.2)
par(mar = marplot)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame20$mop_basic, col = mopame_bas20$col, 
      breaks = mopame_bas20$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame20$mop_simple, col = mopame_sim20$col, 
      breaks = mopame_sim20$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

## 10000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (10,000)", cex = 1.2)
par(mar = marplot)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame10$mop_basic, col = mopame_bas10$col, 
      breaks = mopame_bas10$breaks, add = TRUE)
image(mop_ame10$mop_simple, col = "#000000", add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame10$mop_simple, col = mopame_sim10$col, 
      breaks = mopame_sim10$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

## 5000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (5,000)", cex = 1.2)
par(mar = marplot)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame5$mop_basic, col = mopame_bas5$col, 
      breaks = mopame_bas5$breaks, add = TRUE)
image(mop_ame5$mop_simple, col = "#000000", add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame5$mop_simple, col = mopame_sim5$col, 
      breaks = mopame_sim5$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()

dev.off()
# ------------------------------------------------------------------------------



# Figure S1 ---------------------------------------------------------------------
# colors
eura10 <- plot_raster_inf(mop_ame_eucl10, result = "basic", 
                          color_palette = daright)

eusc10 <- plot_raster_inf(mop_ame_eu_s10, result = "basic", 
                          color_palette = daright)

maha10 <- plot_raster_inf(mop_ame_maha10, result = "basic", 
                          color_palette = daright)

eura5 <- plot_raster_inf(mop_ame_eucl5, result = "basic", 
                         color_palette = daright)

eusc5 <- plot_raster_inf(mop_ame_eu_s5, result = "basic", 
                         color_palette = daright)

maha5 <- plot_raster_inf(mop_ame_maha5, result = "basic", 
                         color_palette = daright)

eura1 <- plot_raster_inf(mop_ame_eucl1, result = "basic", 
                         color_palette = daright)

eusc1 <- plot_raster_inf(mop_ame_eu_s1, result = "basic", 
                         color_palette = daright)

maha1 <- plot_raster_inf(mop_ame_maha1, result = "basic", 
                         color_palette = daright)

# size proportions
cext <- 0.95
cexe <- 0.8
cexl <- 0.7
cexb <- 1.2
ptcex <- 0.5

# initial figure showing example
png("Figures/FigureS1.png", res = 600, width = 166, height = 145, units = "mm")
par(mfrow = c(3, 3), mar = c(4.2, 4.2, 0.65, 0.5), cex = 0.5)

## 10%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = eura10$col[eura10$val_factor], 
     cex = ptcex, xlab = "", ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Euclidean-raw variables)", bty = "n", 
       cex = cext)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = eusc10$col[eusc10$val_factor], 
     cex = ptcex, xlab = "", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = maha10$col[maha10$val_factor],  
     xlab = "", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 10% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

## 5%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = eura5$col[eura5$val_factor], 
     cex = ptcex, xlab = "", ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Euclidean-raw variables)", bty = "n", 
       cex = cext)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = eusc5$col[eusc5$val_factor], 
     cex = ptcex, xlab = "", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = maha5$col[maha5$val_factor],  
     xlab = "", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 5% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

## 1%
## results from euclidean distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = eura1$col[eura1$val_factor], 
     cex = ptcex, xlab = "Minimum temperature", 
     ylab = "Precipitation of wettest month")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Euclidean-raw variables)", bty = "n", 
       cex = cext)
legend(x = -480, y = 850, legend = "Non-analogous", fill = "#010101", bty = "n", 
       cex = cexe)
legend_bar(position = c(-450, 450), col = eura1$col, title = "Dissimilarity", 
           cex = cexb)

## results from euclidean distances and scaled variables
plot(americasv[, c(2, 4)], pch = 16, col = eusc1$col[eusc1$val_factor], 
     cex = ptcex, xlab = "Minimum temperature", ylab = "")
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Euclidean-scaled variables)", bty = "n", 
       cex = cext)

## results from mahalanobis distances and raw variables
plot(americasv[, c(2, 4)], pch = 16, col = maha1$col[maha1$val_factor],  
     xlab = "Minimum temperature", ylab = "", cex = ptcex)
rec <- lapply(block_coor, function(x) {
  rect(x[1], x[2], x[3], x[4], border = colsel)
})
legend("topleft", legend = "MOP 1% (Mahalanobis-raw variables)", bty = "n", 
       cex = cext)

dev.off()
# ------------------------------------------------------------------------------



# Figure S2 --------------------------------------------------------------------
# variables involved
varsn <- names(mop_ame$mop_detailed$towards_low_end)

# color
col0 <- "#740211"

png("Figures/FigureS2.png", res = 600, width = 166, height = 130, units = "mm")
mat <- matrix(1:35, nrow = 5, byrow = TRUE)
layout(mat, widths = c(0.9, rep(10, 6)), heights = c(0.9, rep(10, 4)))

## labes Y
par(mar = martext)
par(cex = 0.7)
plot.new()
for (i in varsn) {
  plot.new(); text(0.5, 0.4, i, cex = 1)
}

## row1
plot.new(); text(0.6, 0.5, "Towards low end", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(boxpamame, col = NA, axes = FALSE)
  image(mop_ame$mop_detailed$towards_low_end[[i]], col = col0, add = TRUE)
  map(add = TRUE, col = colmap)
  box()
}

## row2
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards high end", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(boxpamame, col = NA, axes = FALSE)
  image(mop_ame$mop_detailed$towards_high_end[[i]], col = col0, add = TRUE)
  map(add = TRUE, col = colmap)
  box()
}

## row3
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards low end", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(boxpamecu, col = NA, axes = FALSE)
  image(mop_ecu$mop_detailed$towards_low_end[[i]], col = col0, add = TRUE)
  map(add = TRUE, col = colmap)
  box()
}

## row4
par(mar = martext, cex = 0.7)
plot.new(); text(0.6, 0.5, "Towards high end", cex = 1, srt = 90)

par(mar = marplot, cex = 0.8)
for (i in varsn) {
  plot(boxpamecu, col = NA, axes = FALSE)
  image(mop_ecu$mop_detailed$towards_high_end[[i]], col = col0, add = TRUE)
  map(add = TRUE, col = colmap)
  box()
}

dev.off()
# ------------------------------------------------------------------------------



# Figure S3 --------------------------------------------------------------------
# figure of combination of variables outside ranges
## margins 
marplot <- c(0.1, 0.1, 0.1, 0.1)
martext <- rep(0, 4)

## preparing colors
mopame_low <- plot_raster_inf(mop_amea, result = "towards_low_combined", 
                              color_palette = bluered)

mopame_low5 <- plot_raster_inf(mop_ame5, result = "towards_low_combined", 
                               color_palette = bluered)

mopame_low10 <- plot_raster_inf(mop_ame10, result = "towards_low_combined", 
                                color_palette = bluered)

mopame_low20 <- plot_raster_inf(mop_ame20, result = "towards_low_combined", 
                                color_palette = bluered)

mopame_high <- plot_raster_inf(mop_amea, result = "towards_high_combined", 
                              color_palette = bluered)

mopame_high5 <- plot_raster_inf(mop_ame5, result = "towards_high_combined", 
                               color_palette = bluered)

mopame_high10 <- plot_raster_inf(mop_ame10, result = "towards_high_combined", 
                                color_palette = bluered)

mopame_high20 <- plot_raster_inf(mop_ame20, result = "towards_high_combined", 
                                color_palette = bluered)

## the figure
png("Figures/FigureS3.png", res = 600, width = 160, height = 120, units = "mm")

mat <- matrix(1:15, nrow = 3)
layout(mat, widths = c(1, rep(10, 4)), heights = c(1, 10, 10))

## labels Y
par(mar = martext)
par(cex = 0.6)
plot.new()
plot.new(); text(0.5, 0.5, "Non-analogous towards low end", cex = 1.2, srt = 90)
plot.new(); text(0.5, 0.5, "Non-analogous towards high end", cex = 1.2, srt = 90)

## All data
plot.new(); text(0.5, 0.5, "All data (40,207)", cex = 1.2)

par(mar = marplot)
plot(boxpamame, col = NA, axes = FALSE)
image(mop_amea$mop_detailed$towards_low_combined, col = mopame_low$col, 
      breaks = mopame_low$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_low$legend_text, 
       fill = mopame_low$legend_col, bty = "n", cex = 0.8)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_amea$mop_detailed$towards_high_combined, col = mopame_high$col, 
      breaks = mopame_high$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_high$legend_text, 
       fill = mopame_high$legend_col, bty = "n", cex = 0.8)

## 20000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (20,000)", cex = 1.2)
par(mar = marplot)

par(mar = marplot)
plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame20$mop_detailed$towards_low_combined, col = mopame_low20$col, 
      breaks = mopame_low20$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_low20$legend_text, 
       fill = mopame_low20$legend_col, bty = "n", cex = 0.8)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame20$mop_detailed$towards_high_combined, col = mopame_high20$col, 
      breaks = mopame_high20$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_high20$legend_text, 
       fill = mopame_high20$legend_col, bty = "n", cex = 0.8)

## 10000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (10,000)", cex = 1.2)
par(mar = marplot)

par(mar = marplot)
plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame10$mop_detailed$towards_low_combined, col = mopame_low10$col, 
      breaks = mopame_low10$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_low10$legend_text, 
       fill = mopame_low10$legend_col, bty = "n", cex = 0.8)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame10$mop_detailed$towards_high_combined, col = mopame_high10$col, 
      breaks = mopame_high10$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_high10$legend_text, 
       fill = mopame_high10$legend_col, bty = "n", cex = 0.8)

## 5000
par(mar = martext); plot.new(); text(0.5, 0.5, "Sample (5,000)", cex = 1.2)
par(mar = marplot)

par(mar = marplot)
plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame5$mop_detailed$towards_low_combined, col = mopame_low5$col, 
      breaks = mopame_low5$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_low5$legend_text, 
       fill = mopame_low5$legend_col, bty = "n", cex = 0.8)

plot(boxpamame, col = NA, axes = FALSE)
image(mop_ame5$mop_detailed$towards_high_combined, col = mopame_high5$col, 
      breaks = mopame_high5$breaks, add = TRUE)
map(add = TRUE, col = alpha("gray55", 0.6))
box()
legend("bottomleft", legend = mopame_high5$legend_text, 
       fill = mopame_high5$legend_col, bty = "n", cex = 0.8)

dev.off()
# ------------------------------------------------------------------------------
