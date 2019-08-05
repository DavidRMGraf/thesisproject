# script diffusion map complete - 24.07.19
#with otus, physical data (and optionally nutrients!)

### this part was taken from the universal header! (START) ----
rm(list = ls())
graphics.off()

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
phys_oce <- readRDS("physical_oceanography_data_all_years.rds")

#remove autocorrelated values from the physical oceanography dataset:
phys_oce.sub.all <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                "NO3_mumol_l", "NO2_mumol_l", "SiOH4_mumol_l", "PO4_mumol_l"))
phys_oce.sub.phy <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit"))

# cases need to be complete.cases AND the duplicates need to be excluded:
columns2keep.all <- complete.cases(phys_oce.sub.all) & phys_oce$keep == 1   
columns2keep.phy <- complete.cases(phys_oce.sub.phy) & phys_oce$keep == 1
# reduce
phys_oce.sub.all <- phys_oce.sub.all[columns2keep.all,]
phys_oce.sub.phy <- phys_oce.sub.phy[columns2keep.phy,]

# get sequences
sequ <- readRDS("sequ_all.rds")

# apply 0.05 percent threshold
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")

dim(sequ)

sequ.sub.all <- sequ[columns2keep.all,]
sequ.sub.phy <- sequ[columns2keep.phy,]

# get station names from the phys_oce datasheet:
stat_names <- subset(phys_oce, select = c("Proben_ID_intern", "date", "depth",
                                          "year", "latitude", "longitude"))
# reduce
stat_names.all <- stat_names[columns2keep.all,]
stat_names.phy <- stat_names[columns2keep.phy,]
### this part was taken from the universal header! (END) ----
# remove unwanted filters, original datasets 
rm(columns2keep.all, columns2keep.phy, sequ, stat_names, phys_oce)
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
data.all <- similarity(as.matrix(f.clr.all))
data.phy <- similarity(as.matrix(f.clr.phy))

data.all <- simil_reducer(data.all)
data.phy <- simil_reducer(data.phy)


lap.all <- matrixLaplacian::matrixLaplacian(data.all, plot2D = F, plot3D = F)
lap.phy <- matrixLaplacian::matrixLaplacian(data.phy, plot2D = F, plot3D = F)

lap_mat.all <- lap.all$LaplacianMatrix
lap_mat.phy <- lap.phy$LaplacianMatrix

elm.all <- eigen(lap_mat.all)
elm.phy <- eigen(lap_mat.phy)

# construct data.frames from eigenvectors with lowest values from each elm result:

low.all <- data.frame("minus_1" = elm.all$vectors[,ncol(data.all)-1], "minus_2" = elm.all$vectors[,ncol(data.all)-2],
                      "minus_3" = elm.all$vectors[,ncol(data.all)-3], "minus_4" = elm.all$vectors[,ncol(data.all)-4],
                      "minus_5" = elm.all$vectors[,ncol(data.all)-5], "minus_6" = elm.all$vectors[,ncol(data.all)-6],
                      "minus_7" = elm.all$vectors[,ncol(data.all)-7], "minus_8" = elm.all$vectors[,ncol(data.all)-8],
                      year = as.factor(stat_names.all$year),
                      depth = -as.numeric(stat_names.all$depth),
                      count = 1:ncol(data.all))
low.phy <- data.frame("minus_1" = elm.phy$vectors[,ncol(data.phy)-1], "minus_2" = elm.phy$vectors[,ncol(data.phy)-2],
                      "minus_3" = elm.phy$vectors[,ncol(data.phy)-3], "minus_4" = elm.phy$vectors[,ncol(data.phy)-4],
                      "minus_5" = elm.phy$vectors[,ncol(data.phy)-5], "minus_6" = elm.phy$vectors[,ncol(data.phy)-6],
                      "minus_7" = elm.phy$vectors[,ncol(data.phy)-7], "minus_8" = elm.phy$vectors[,ncol(data.phy)-8],
                      year = as.factor(stat_names.phy$year),
                      depth = -as.numeric(stat_names.phy$depth),
                      count = 1:ncol(data.phy))

ggplot(low.all, aes(x = minus_1 , y = minus_2, col = depth))+
  geom_point()+
  labs(title = "physical+nutrient subset",
       subtitle = "seems to pick up time signal!")

ggplot(low.phy, aes(x = minus_1 , y = minus_2, col = depth))+
  geom_point()+
  labs(title = "physical subset",
       subtitle = "seems to pick up time signal!")



ggplot(low.all, aes(x = minus_2 , y = minus_3, col =depth))+
  geom_point()+
  labs(subtitle = "seems to pick up time signal!")



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
