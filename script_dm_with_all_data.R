# script diffusion map complete - 24.07.19
#with otus, physical data (and optionally nutrients!)

rm(list = ls())
graphics.off()

## packages ---------------------------------------------------

## functions --------------------------------------------------
similarity <- function(input_data){
  if(!is.matrix(input_data)) stop('input_data must be a matrix')
  simil <- matrix(data = NA, nrow = nrow(input_data), ncol = nrow(input_data))
  for (i in 1:nrow(input_data)){
    for (j in 1:nrow(input_data)){
      simil[i,j] <- dist(rbind(input_data[i, ], input_data[j, ]), method = "euclidean")
    }
  }
  simil <- 1/simil
  diag(simil) <- 0
  return(simil)
}
simil_reducer <- function(simil){
  if(!is.matrix(simil)) stop('Input must be a matrix')
  simil_red <- matrix(data = NA, nrow = nrow(simil), ncol = nrow(simil))
  for (i in 1:nrow(simil)){
    for (j in 1:nrow(simil)){
      if (simil[i, j] < min(c(sort(simil[i, ], decreasing = T)[10], sort(simil[, j], decreasing = T)[10]))){
        simil_red[i, j] <- 0
      }else{
        simil_red[i, j] <- simil[i, j]
      }
    }
  }
  return(simil_red)
}

## data -------------------------------------------------------
# get physical oceanography data
phys_oce <- readRDS("physical_oceanography_data_all_years.rds")

#remove autocorrelated values from the physical oceanography dataset:
phys_oce.sub.all <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                "NO3_mumol_l", "NO2_mumol_l", "SiOH4_mumol_l", "PO4_mumol_l"))
phys_oce.sub.phy <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit"))

## cases need to be complete.cases AND the duplicates need to be excluded:
columns2keep.all <- complete.cases(phys_oce.sub.all) & phys_oce$keep == 1   
columns2keep.phy <- complete.cases(phys_oce.sub.phy) & phys_oce$keep == 1

phys_oce.sub.all <- phys_oce.sub.all[columns2keep.all,]
phys_oce.sub.phy <- phys_oce.sub.phy[columns2keep.phy,]

# get sequences, 0.05 percent threshold applied
sequ <- readRDS("sequences_thresh_applied.rds")
dim(sequ)

sequ.sub.all <- sequ[columns2keep.all,]
sequ.sub.phy <- sequ[columns2keep.phy,]

# get station names from the phys_oce datasheet:
stat_names <- subset(phys_oce, select = c("Proben_ID_intern", "date", "depth",
                                          "year", "latitude", "longitude"))

stat_names.all <- stat_names[columns2keep.all,]
stat_names.phy <- stat_names[columns2keep.phy,]


## DM --------------------------------------------------------
# variance-stabilizing transformation:

f.n0.all <- zCompositions::cmultRepl(sequ.sub.all, method="CZM", label = 0)
f.n0.phy <- zCompositions::cmultRepl(sequ.sub.phy, method="CZM", label = 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# physical data oder nutrients+physical data dazu?
# f.n0.all <- cbind(f.n0.all, phys_oce.sub.all)
# f.n0.phy <- cbind(f.n0.phy, phys_oce.sub.phy)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f.clr.all <- CoDaSeq::codaSeq.clr(f.n0.all, samples.by.row = T)
f.clr.phy <- CoDaSeq::codaSeq.clr(f.n0.phy, samples.by.row = T)

# DM (Thilo's method)
data.all <- similarity(as.matrix(f.n0.all))
data.phy <- similarity(as.matrix(f.n0.phy))

data.all <- simil_reducer(data.all)
data.phy <- simil_reducer(data.phy)


lap.all <- matrixLaplacian::matrixLaplacian(data.all, plot2D = F, plot3D = F)
lap.phy <- matrixLaplacian::matrixLaplacian(data.phy, plot2D = F, plot3D = F)

lap_mat.all <- lap.all$LaplacianMatrix
lap_mat.phy <- lap.phy$LaplacianMatrix

elm.all <- eigen(lap_mat.all)
elm.phy <- eigen(lap_mat.phy)

## map plot --------------------------------------------------
library(rgdal)                                                                                                      
library(raster)
library(ggplot2)
library(rworldxtra)
library(rgeos)

stat_coord.all <- stat_names.all[,c("latitude","longitude")] # select longitude and latitude columns
#stat_coord <- na.omit(stat_coord.all) # remove missing values
coordinates(stat_coord.all) <- ~longitude+latitude # spatial coordinate information 


# get world map
wm <- rworldmap::getMap(resolution ="high")
wm <- crop(wm, extent(-20, 20, 76, 81)) # adjust extent of map to sampling region# simple plot

plot(wm, axes=T) 
points(coordinates(stat_coord.all)) 
