## script maps in R

rm(list=ls())
graphics.off()
# packages ---------------------------
library(raster)
library(readxl)
library(ggplot2)
# map plotting -----------------------
all_oce <- shapefile("ne_10m_ocean.shp")
stat_dat <- read_xlsx("Sequences_Hausgarten_station_data_revised.xlsx")

# Basic plot of this shape file:

plot(all_oce,
     xlim = c(-20, 20),
     ylim = c(77, 77.5),
     col = "white",
     bg = "black")

stat_dat <- stat_dat[order(stat_dat$year),]
points(stat_dat$longitude, stat_dat$latitude,
       col = stat_dat$year-2008,
       pch = 20)
legend("bottom", legend = as.character(unique(stat_dat$year)), pch = 20, col = unique(stat_dat$year)-2008, horiz = T)
