# script to calculate the physical distance matrix from the supplementary data sheet
rm(list=ls())
graphics.off()
## packages ----------------
library(geosphere)

## script ------------------
stations <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

stations <- stations[stations$latitude<360 & !is.na(stations$latitude),]
# uni <- unique(stations$station)
# statcoord <- stations[stations$station == uni[1], ][1,]
# 
# for (i in 2:length(uni)){
#   statcoord <- rbind(statcoord, stations[stations$station == uni[i], ][1,])
# }
# statcoord <- statcoord[, 7:9]
# 
# mat <- distm(statcoord[, 2:3], statcoord[, 2:3], fun = distVincentyEllipsoid)
# rownames(mat) <- statcoord[,1]

mat <- distm(stations[, 8:9], stations[, 8:9], fun = distVincentyEllipsoid)
