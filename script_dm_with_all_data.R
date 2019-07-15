# script diffusion map complete with otus, physical data (and optionally nutrients!)
rm(list = ls())
graphics.off()

## packages ---------------------------------------------------
library(readxl)

## functions --------------------------------------------------
threshapply <- function(input_data, method){
  if(!is.character(method)) stop('method must be character')
  if(!is.matrix(input_data)) stop('input_data must be a matrix')
  if (method=="90 percent"){
    for (i in 1:nrow(input_data)){
      output_data <- input_data
      output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.90]] <- 0
    }
  }else if(method == "95 percent"){
    for (i in 1:nrow(input_data)){
      output_data <- input_data
      output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.95]] <- 0
    }
  }else if(method == "99 percent"){
    for (i in 1:nrow(input_data)){
      output_data <- input_data
      output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.99]] <- 0
    }
  }else if(method == "0.05 percent"){
    keep.cols <- colSums(input_data)/sum(input_data)>=5e-04
    output_data <- input_data[, keep.cols]
  }else if(method == "0.005 percent"){
    keep.cols <- colSums(input_data)/sum(input_data)>=5e-05
    output_data <- input_data[, keep.cols]
  }
  output_data <- output_data[, colSums(output_data)!=0]
  return(output_data)
}
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
sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";", header = F)
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")

phys <- read.csv("physical_data_matrix.csv", header = T)[, 12:22]

for(i in 1:ncol(phys)){print(colnames(phys)[i]); print(sum(is.na(phys[,i])))}


data <- as.matrix(cbind(sequ, phys))

stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

## DM --------------------------------------------------------
# variance-stabilizing transformation:
f.n0 <- zCompositions::cmultRepl(sequ, method="CZM", label = 0)
f.clr <- CoDaSeq::codaSeq.clr(f.n0, samples.by.row = T)

# DM (Thilo's method)
data <- similarity(as.matrix(f.n0))
data <- simil_reducer(data)

lap <- matrixLaplacian(data, plot2D = F, plot3D = F)
lap_mat <- lap$LaplacianMatrix

elm <- eigen(lap_mat)

## map plot --------------------------------------------------
library(rgdal)                                                                                                      
library(raster)
library(ggplot2)
library(rworldxtra)
library(rgeos)

stat_coord <- stat_names[,c("latitude","longitude")] # select longitude and latitude columns
stat_coord <- na.omit(stat_coord) # remove missing values
coordinates(stat_coord) <- ~longitude+latitude # spatial coordinate information 


# get world map
wm <- rworldmap::getMap(resolution ="high")
wm <- crop(wm, extent(-20, 20, 76, 81)) # adjust extent of map to sampling region# simple plot

plot(wm, axes=T) 
points(coordinates(stat_coord), pch=stat_names$year-2008, col="red") 
