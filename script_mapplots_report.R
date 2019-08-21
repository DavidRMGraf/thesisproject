# maps for the thesis report
## resource:
# https://www.molecularecologist.com/2012/09/making-maps-with-r/
rm(list=ls())
graphics.off()

setwd("~/Studium/19SS/BA/thesisproject")
library(maps)
library(mapdata)

library(Cairo)
# get analysis data:
source("script_analysis_all_data.r")
rm(list = ls()[-(38:40)])

setwd("~/Studium/19SS/BA/ba_thesis_report/map_plots")

mapplot <- function(dat.in){
  map("worldHires","Norway:Svalbard", xlim = c(-20,15), ylim = c(76, 81), col="gray90", fill=TRUE)
  map("worldHires","Greenland", xlim = c(-20,15), ylim = c(76, 81), col="gray90", fill=TRUE, add = TRUE)
  points(dat.in$longitude, dat.in$latitude, pch=19, col="red", cex=0.3)
  rect(xleft = -5, ybottom = 78.5, xright = 11, ytop = 80, col = NA, border = "black", lty = "dotted", lwd = 1)
  map.axes()
}

mapplot(stat_names.2014)

## 2014 map: -------
Cairo(file="subset_2014.png", 
      type="png",
      units="in", 
      width=5*2, 
      height=4*2, 
      pointsize=12*2, 
      dpi=72)
mapplot(stat_names.2014)
dev.off()

## 2016 map: -------
Cairo(file="subset_2016.png", 
      type="png",
      units="in", 
      width=5*2, 
      height=4*2, 
      pointsize=12*2, 
      dpi=72)
mapplot(stat_names.2016)
dev.off()

## long map: -------
Cairo(file="subset_long.png", 
      type="png",
      units="in", 
      width=5*2, 
      height=4*2, 
      pointsize=12*2, 
      dpi=72)
mapplot(stat_names.long)
dev.off()

