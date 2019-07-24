## script to perform the regressions on distance matrices - 24.07.19
rm(list=ls())
graphics.off()

library(ecodist)
library(geosphere)
library(vegan)

#script -----------------------------------------------------------
# get physical oceanography data
phys_oce <- readRDS("physical_oceanography_data_all_years.rds")

# get sequences, 0.05 percent threshold applied
sequ <- readRDS("sequences_thresh_applied.rds")

# determine which columns to keep, based on physical dataset:
col2keep <- complete.cases(phys_oce)

# reducing physical oceanography dataset:
phys_oce <- phys_oce[col2keep,]
range(phys_oce$year)
unique(phys_oce$year)

# reducing sequence dataset:
sequ <- sequ[col2keep,]

## physical distance matrix: --------------------------------------
X <- phys_oce$longitude
Y <- phys_oce$latitude
points <- as.matrix(data.frame(X,Y))

phys.dist.m <- geosphere::distm(points, fun = distVincentyEllipsoid)

## temporal distance matrix: --------------------------------------
temp.m  <- as.matrix(dist(phys_oce$year, method="euclidean", diag = FALSE, upper = FALSE))  

## bray-curtis similarity matrix on sequences: --------------------
bray.dis <- vegan::vegdist(sequ, method='bray')
bray.sim <- 1-as.matrix(bray.dis)

## MLM: -----------------------------------------------------------
d <- list(bray.sim, phys.dist.m, temp.m)
names(d) <- c('bray.sim','space','time')
lower  <- lapply(d, ecodist::lower)
d      <- as.data.frame(lapply(lower, c))

m <- ecodist::MRM(bray.sim ~ space + time, data = d); m

