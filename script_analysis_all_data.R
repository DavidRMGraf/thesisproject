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
library(FactoMineR)
library(factoextra)

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

# sum up nitrogen species:
phys_oce$nitrogen_species <- phys_oce$NO2_mumol_l + phys_oce$NO3_mumol_l

#remove autocorrelated values from the physical oceanography dataset:
phys_oce.sub.2014 <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                "nitrogen_species", "SiOH4_mumol_l", "PO4_mumol_l"))
phys_oce.sub.long <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                 "nitrogen_species", "SiOH4_mumol_l", "PO4_mumol_l"))
phys_oce.sub.2016 <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit", "icecover"))

# cases need to be complete.cases AND the duplicates need to be excluded:
# only values from 2014 in the whole Fram Strait:
columns2keep.2014 <- complete.cases(phys_oce.sub.2014) & phys_oce$keep == 1 & phys_oce$year == 2014
# all years except 2016, shallower than 30.1 m, in HAUSGARTEN:
columns2keep.long <- complete.cases(phys_oce.sub.long) & phys_oce$keep == 1 & phys_oce$depth <= 30 & phys_oce$latitude >= 78.5 & phys_oce$latitude <=80 & phys_oce$longitude >= -5 & phys_oce$longitude <= 11
# # borders of HG after: https://www.awi.de/en/science/biosciences/deep-sea-ecology-and-technology/observatories/lter-observatory-hausgarten.html
# everything available from 2016:
columns2keep.2016 <- complete.cases(phys_oce.sub.2016) & phys_oce$keep == 1 & phys_oce$year == 2016

# reduce
phys_oce.sub.2014 <- phys_oce.sub.2014[columns2keep.2014,]
phys_oce.sub.long <- phys_oce.sub.long[columns2keep.long,]
phys_oce.sub.2016 <- phys_oce.sub.2016[columns2keep.2016,]

# get sequences
sequ <- readRDS("sequ_all.rds")
sequ <- t(sequ)
# apply 0.05 percent threshold
dim(sequ)

sequ.sub.2014 <- threshapply(sequ[columns2keep.2014,], "0.05 percent")
sequ.sub.long <- threshapply(sequ[columns2keep.long,], "0.05 percent")
sequ.sub.2016 <- threshapply(sequ[columns2keep.2016,], "0.05 percent")
# sequ.sub.all <- sequ[columns2keep.all,]
# sequ.sub.phy <- sequ[columns2keep.phy,]
# sequ.sub.HG <- sequ[columns2keep.HG,]

# get station names from the phys_oce datasheet:
stat_names <- subset(phys_oce, select = c("Proben_ID_intern", "date", "depth",
                                          "year", "latitude", "longitude", "icecover"))

# reduce
stat_names.2014 <- stat_names[columns2keep.2014,]
stat_names.long <- stat_names[columns2keep.long,]
stat_names.2016 <- stat_names[columns2keep.2016,]
# stat_names.all <- stat_names[columns2keep.all,]
# stat_names.phy <- stat_names[columns2keep.phy,]
# stat_names.HG <- stat_names[columns2keep.HG,]

### this part was taken from the universal header! (END) ----
# remove unwanted filters, original datasets 
rm(columns2keep.2014, columns2keep.2016, columns2keep.long,
   sequ, stat_names, phys_oce,
   threshapply)

## triggers:--------------------------------------------------

# calculate correlations?
calc.correl <- FALSE
# plot worldmap?
plot.wm <- FALSE

## Transformation -------------------------------------------
# remove zeros in sequence subset:
f.n0.2014 <- zCompositions::cmultRepl(sequ.sub.2014, method="CZM", label = 0)
f.n0.long <- zCompositions::cmultRepl(sequ.sub.long, method="CZM", label = 0)
f.n0.2016 <- zCompositions::cmultRepl(sequ.sub.2016, method="CZM", label = 0)

# variance-stabilizing transformation:
f.clr.2014 <- CoDaSeq::codaSeq.clr(f.n0.2014, samples.by.row = T)
f.clr.long <- CoDaSeq::codaSeq.clr(f.n0.long, samples.by.row = T)
f.clr.2016 <- CoDaSeq::codaSeq.clr(f.n0.2016, samples.by.row = T)

## PCA ------------------------------------------------------
pca.2014 <- PCA(f.clr.2014, scale.unit = F, ncp = 5, graph = FALSE)
pca.long <- PCA(f.clr.long, scale.unit = F, ncp = 5, graph = FALSE)
pca.2016 <- PCA(f.clr.2016, scale.unit = F, ncp = 5, graph = FALSE)

# get coordinates in the PCA for each station:
coords.2014 <- as.data.frame.matrix(pca.2014$ind$coord)
coords.long <- as.data.frame.matrix(pca.long$ind$coord)
coords.2016 <- as.data.frame.matrix(pca.2016$ind$coord)

## DM -------------------------------------------------------
data.2014 <- similarity(as.matrix(f.clr.2014))
data.long <- similarity(as.matrix(f.clr.long))
data.2016 <- similarity(as.matrix(f.clr.2016))

data.2014 <- simil_reducer(data.2014)
data.long <- simil_reducer(data.long)
data.2016 <- simil_reducer(data.2016)

lap.2014 <-  matrixLaplacian(data.2014, plot2D = F, plot3D = F)
lap.long <-  matrixLaplacian(data.long, plot2D = F, plot3D = F)
lap.2016 <-  matrixLaplacian(data.2016, plot2D = F, plot3D = F)

lap_mat.2014 <- lap.2014$LaplacianMatrix
lap_mat.long <- lap.long$LaplacianMatrix
lap_mat.2016 <- lap.2016$LaplacianMatrix

elm.2014 <- eigen(lap_mat.2014)
elm.long <- eigen(lap_mat.long)
elm.2016 <- eigen(lap_mat.2016)

# construct data.frames from eigenvectors with lowest values from each elm result:
low.2014 <- data.frame("minus_1" = elm.2014$vectors[,ncol(data.2014)-1], "minus_2" = elm.2014$vectors[,ncol(data.2014)-2],
                       "minus_3" = elm.2014$vectors[,ncol(data.2014)-3], "minus_4" = elm.2014$vectors[,ncol(data.2014)-4],
                       "minus_5" = elm.2014$vectors[,ncol(data.2014)-5], "minus_6" = elm.2014$vectors[,ncol(data.2014)-6],
                       "minus_7" = elm.2014$vectors[,ncol(data.2014)-7], "minus_8" = elm.2014$vectors[,ncol(data.2014)-8])
low.long <- data.frame("minus_1" = elm.long$vectors[,ncol(data.long)-1], "minus_2" = elm.long$vectors[,ncol(data.long)-2],
                       "minus_3" = elm.long$vectors[,ncol(data.long)-3], "minus_4" = elm.long$vectors[,ncol(data.long)-4],
                       "minus_5" = elm.long$vectors[,ncol(data.long)-5], "minus_6" = elm.long$vectors[,ncol(data.long)-6],
                       "minus_7" = elm.long$vectors[,ncol(data.long)-7], "minus_8" = elm.long$vectors[,ncol(data.long)-8])
low.2016 <- data.frame("minus_1" = elm.2016$vectors[,ncol(data.2016)-1], "minus_2" = elm.2016$vectors[,ncol(data.2016)-2],
                       "minus_3" = elm.2016$vectors[,ncol(data.2016)-3], "minus_4" = elm.2016$vectors[,ncol(data.2016)-4],
                       "minus_5" = elm.2016$vectors[,ncol(data.2016)-5], "minus_6" = elm.2016$vectors[,ncol(data.2016)-6],
                       "minus_7" = elm.2016$vectors[,ncol(data.2016)-7], "minus_8" = elm.2016$vectors[,ncol(data.2016)-8])
## correlations -------
#' leaving out the first (smallest non-zero eigenvalue) and the corresponding eigenvector,
#' each eigenvector is correlated with each column in the original dataset to find correlations
#' between the eigenvectors and the predictors, determining what feature the eigenvector picks up

# here, f.clr should be replaced with the untransformed values f.n0 as Katja confirmed!
if (calc.correl){
  ## correlations 2014:
  cor.coef.2014 <- data.frame("eig.vec.ind" = sort(rep((length(elm.2014$values)-1):1, 10), decreasing = T),
                             "order.ev" = rep(1:10, length(elm.2014$values)-1), "p-value" = NA,
                             "best.predictor" = NA, "cor.coeff" = NA)
  coeffic <- NA
  pval <- NA
  for(i in (length(elm.2014$values)-1):1){
    for(j in 1:dim(f.n0.input.2014)[2]){
      coeffic[j] <- cor.test(elm.2014$vectors[, i], f.n0.input.2014[,j])[[4]]
      pval[j] <- cor.test(elm.2014$vectors[, i], f.n0.input.2014[,j])[[3]]
    }
    logi <- order(abs(coeffic), decreasing = T)[1:10]
    cor.coef.2014$best.predictor[cor.coef.2014$eig.vec.ind == i] <- paste("Var. ", colnames(f.n0.input.2014)[logi])
    cor.coef.2014$cor.coeff[cor.coef.2014$eig.vec.ind == i] <- coeffic[logi]
    cor.coef.2014$p.value[cor.coef.2014$eig.vec.ind == i] <- pval[logi]
  }
  
  ## correlations long:
  cor.coef.long <- data.frame("eig.vec.ind" = sort(rep((length(elm.long$values)-1):1, 10), decreasing = T),
                              "order.ev" = rep(1:10, length(elm.long$values)-1), "p-value" = NA,
                              "best.predictor" = NA, "cor.coeff" = NA)
  coeffic <- NA
  pval <- NA
  for(i in (length(elm.long$values)-1):1){
    for(j in 1:dim(f.n0.input.long)[2]){
      coeffic[j] <- cor.test(elm.long$vectors[, i], f.n0.input.long[,j])[[4]]
      pval[j] <- cor.test(elm.long$vectors[, i], f.n0.input.long[,j])[[3]]
    }
    logi <- order(abs(coeffic), decreasing = T)[1:10]
    cor.coef.long$best.predictor[cor.coef.long$eig.vec.ind == i] <- paste("Var. ", colnames(f.n0.input.long)[logi])
    cor.coef.long$cor.coeff[cor.coef.long$eig.vec.ind == i] <- coeffic[logi]
    cor.coef.long$p.value[cor.coef.long$eig.vec.ind == i] <- pval[logi]
  }
  
  ## correlations 2016:
  cor.coef.2016 <- data.frame("eig.vec.ind" = sort(rep((length(elm.2016$values)-1):1, 10), decreasing = T),
                              "order.ev" = rep(1:10, length(elm.2016$values)-1), "p-value" = NA,
                              "best.predictor" = NA, "cor.coeff" = NA)
  coeffic <- NA
  pval <- NA
  for(i in (length(elm.2016$values)-1):1){
    for(j in 1:dim(f.n0.input.2016)[2]){
      coeffic[j] <- cor.test(elm.2016$vectors[, i], f.n0.input.2016[,j])[[4]]
      pval[j] <- cor.test(elm.2016$vectors[, i], f.n0.input.2016[,j])[[3]]
    }
    logi <- order(abs(coeffic), decreasing = T)[1:10]
    cor.coef.2016$best.predictor[cor.coef.2016$eig.vec.ind == i] <- paste("Var. ", colnames(f.n0.input.2016)[logi])
    cor.coef.2016$cor.coeff[cor.coef.2016$eig.vec.ind == i] <- coeffic[logi]
    cor.coef.2016$p.value[cor.coef.2016$eig.vec.ind == i] <- pval[logi]
  }
  
  # ## correlations dummy:
  # cor.coef.dummy <- data.frame("eig.vec.ind" = sort(rep((length(elm.dummy$values)-1):1, 10), decreasing = T),
  #                             "order.ev" = rep(1:10, length(elm.dummy$values)-1), "p-value" = NA,
  #                             "best.predictor" = NA, "cor.coeff" = NA)
  # coeffic <- NA
  # pval <- NA
  # for(i in (length(elm.dummy$values)-1):1){
  #   for(j in 1:dim(f.n0.dummy)[2]){
  #     coeffic[j] <- cor.test(elm.dummy$vectors[, i], f.n0.dummy[,j])[[4]]
  #     pval[j] <- cor.test(elm.dummy$vectors[, i], f.n0.dummy[,j])[[3]]
  #   }
  #   logi <- order(abs(coeffic), decreasing = T)[1:10]
  #   cor.coef.dummy$best.predictor[cor.coef.dummy$eig.vec.ind == i] <- paste("Var. ", colnames(f.n0.dummy)[logi])
  #   cor.coef.dummy$cor.coeff[cor.coef.dummy$eig.vec.ind == i] <- coeffic[logi]
  #   cor.coef.dummy$p.value[cor.coef.dummy$eig.vec.ind == i] <- pval[logi]
  # }
  
  cor.coef.2014[cor.coef.2014$order.ev == 1,]
}

## map plot --------------------------------------------------
if (plot.wm){
  library(rgdal)                                                                                                      
  library(raster)
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
}