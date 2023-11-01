# function to make raster and scatter plots share color palette and scale
raster_shared_colors <- function(raster_layers, color_palette = gray.colors, 
                                 raster_plot = TRUE, rev = TRUE) {
  nl <- terra::nlyr(raster_layers)
  
  dat_euc <- terra::as.data.frame(raster_layers)
  
  if (raster_plot) {
    dat_euc1 <- lapply(dat_euc, function(x) {unique(sort(x))})
    nps <- sapply(dat_euc1, length)
    dat_euc1 <- as.factor(unlist(dat_euc1))
  } else {
    dat_euc <- terra::as.data.frame(raster_layers)
    dat_euc1 <- as.factor(unlist(dat_euc))
    np <- nrow(dat_euc)
  }
  
  
  if (rev) {
    col_euc <- rev(color_palette(length(levels(dat_euc1))))
  } else {
    col_euc <- color_palette(length(levels(dat_euc1)))
  }
  
  all_col <- col_euc[dat_euc1]
  
  if (raster_plot) {
    nps <- sapply(1:length(nps), function(x) sum(nps[1:x]))
    seqcol <- c(0, nps)
  } else {
    seqcol <- seq(from = 0, to = length(all_col), length.out = (nl + 1))
  }
  
  lapply(1:nl, function(x) {
    all_col[(seqcol[x] + 1):seqcol[(x + 1)]]
  })
}


# boxplot blocks
block_bxplot <- function(list_blocks, column = 1, ylim = NULL, axis1 = TRUE, 
                         axis2 = TRUE, colbox = "gray75") {
  
  beuc10 <- lapply(list_blocks, function(x) {x[, column]})
  beuc10 <- data.frame(lapply(beuc10, "length<-", max(lengths(beuc10))))
  colnames(beuc10) <- paste("Block", 1:6)
  
  xlab <- ifelse(axis1, "Distance", "")
  
  bx <- boxplot(beuc10, ylab = "", horizontal = TRUE, axes = FALSE,
                xlab = xlab, outline = FALSE, lwd = 0.5, col = colbox, 
                pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5),
                ylim = ylim)
  if (axis1) {
    axis(1, labels = TRUE, lwd = 0.5)
    axis(2, labels = FALSE, lwd = 0.5)
  } 
  
  if (axis2) {
    axis(1, labels = FALSE, lwd = 0.5)
    axis(2, labels = bx$names, at = 1:6, lwd = 0.5, las = 2)
  }
  
  if (!axis1 & !axis2) {
    axis(1, labels = FALSE, lwd = 0.5)
    axis(2, labels = FALSE, lwd = 0.5)
  }
  
  box(lwd = 0.8)
}


# plot results and blocks
block_scatter <- function(xy, block_coor, ptcex = 1, axis1 = TRUE, axis2 = TRUE, 
                          ptcol = "gray75", colsel = "red") {
  
  xlab <- ifelse(axis1, "Minimum temperature", "")
  ylab <- ifelse(axis2, "Precipitation of driest month", "")
  
  plo <- plot(xy, pch = 16, col = ptcol, cex = ptcex, axes = FALSE,
              xlab = xlab, ylab = ylab)
  rec <- lapply(block_coor, function(x) {
    rect(x[1], x[2], x[3], x[4], border = colsel)
  })
  
  if (axis1) {
    axis(1, labels = TRUE, lwd = 0.5)
    axis(2, labels = FALSE, lwd = 0.5)
  } 
  
  if (axis2) {
    axis(1, labels = FALSE, lwd = 0.5)
    axis(2, labels = TRUE, lwd = 0.5, las = 2)
  }
  
  if (!axis1 & !axis2) {
    axis(1, labels = FALSE, lwd = 0.5)
    axis(2, labels = FALSE, lwd = 0.5)
  }
  
  box(lwd = 0.8)
}
