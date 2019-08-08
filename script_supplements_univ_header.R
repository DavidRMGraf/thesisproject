# script diffusion map complete - 24.07.19
#with otus, physical data (and optionally nutrients!)

### this part was taken from the universal header! (START) ----
# MODIFIED
rm(list = ls())
graphics.off()

## packages
library(readxl)
library(zCompositions) #' to remove non-zero entries from the raw data matrix
library(BiocParallel)
library(curl)
library(CoDaSeq) #' installed from tarball from ggloor's github repo on CoDaSeq
library(robCompositions) #' to calculate the Aitchison distance matrix
library(ggbiplot)
library(matrixLaplacian)
library(destiny)
library(rgl)
library(ggplot2)

## functions
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

## data
# get physical oceanography data
phys_oce <- readRDS("physical_oceanography_data_with_ice.rds")

# add 2 to the temperature data to make it compatible:
phys_oce$temp_deg_c <- phys_oce$temp_deg_c+2

#remove autocorrelated values from the physical oceanography dataset:
phys_oce.sub.all <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                "NO3_mumol_l", "NO2_mumol_l", "SiOH4_mumol_l", "PO4_mumol_l", "icecover"))
phys_oce.sub.phy <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit", "icecover"))

# cases need to be complete.cases AND the duplicates need to be excluded:
columns2keep.all <- complete.cases(phys_oce.sub.all) & phys_oce$keep == 1
columns2keep.phy <- complete.cases(phys_oce.sub.phy) & phys_oce$keep == 1 

columns2keep.2014 <- columns2keep.all & phys_oce$year == 2014
columns2keep.long <- columns2keep.all & phys_oce$depth <= 30
columns2keep.2016 <- columns2keep.phy & phys_oce$year == 2016

# columns2keep.HG <- columns2keep.all & phys_oce$latitude >= 78.5 & phys_oce$latitude <=80 & phys_oce$longitude >= -5 & phys_oce$longitude <= 11
# # borders of HG after: https://www.awi.de/en/science/biosciences/deep-sea-ecology-and-technology/observatories/lter-observatory-hausgarten.html

# reduce
phys_oce.sub.2014 <- phys_oce.sub.all[columns2keep.2014,]
phys_oce.sub.long <- phys_oce.sub.all[columns2keep.long,]
phys_oce.sub.2016 <- phys_oce.sub.phy[columns2keep.2016,]

# phys_oce.sub.HG <- phys_oce.sub.all[columns2keep.HG,]
# phys_oce.sub.all <- phys_oce.sub.all[columns2keep.all,]
# phys_oce.sub.phy <- phys_oce.sub.phy[columns2keep.phy,]

# get sequences
sequ <- readRDS("sequ_all.rds")
sequ <- t(sequ)

dim(sequ)

sequ.sub.2014 <- sequ[columns2keep.2014,]
sequ.sub.long <- sequ[columns2keep.long,]
sequ.sub.2016 <- sequ[columns2keep.2016,]
# sequ.sub.all <- sequ[columns2keep.all,]
# sequ.sub.phy <- sequ[columns2keep.phy,]
# sequ.sub.HG <- sequ[columns2keep.HG,]

# apply 0.05 percent threshold
sequ.sub.2014 <- threshapply(sequ.sub.2014, "0.05 percent")
sequ.sub.long <- threshapply(sequ.sub.long, "0.05 percent")
sequ.sub.2016 <- threshapply(sequ.sub.2016, "0.05 percent")
# sequ.sub.all <- threshapply(sequ.sub.all, "0.05 percent")
# sequ.sub.phy <- threshapply(sequ.sub.phy, "0.05 percent")
# sequ.sub.HG <- threshapply(sequ.sub.HG, "0.05 percent")



# get station names from the phys_oce datasheet:
stat_names <- subset(phys_oce, select = c("Proben_ID_intern", "date", "depth",
                                          "year", "latitude", "longitude"))
# reduce
stat_names.2014 <- stat_names[columns2keep.2014,]
stat_names.long <- stat_names[columns2keep.long,]
stat_names.2016 <- stat_names[columns2keep.2016,]
# stat_names.all <- stat_names[columns2keep.all,]
# stat_names.phy <- stat_names[columns2keep.phy,]
# stat_names.HG <- stat_names[columns2keep.HG,]

### this part was taken from the universal header! (END) ----
