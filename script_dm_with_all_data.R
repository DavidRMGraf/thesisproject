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
phys_oce <- readRDS("physical_oceanography_data_all_years.rds")

# add 2 to the temperature data to make it compatible:
phys_oce$temp_deg_c <- phys_oce$temp_deg_c+2

#remove autocorrelated values from the physical oceanography dataset:
phys_oce.sub.all <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit",
                                                "NO3_mumol_l", "NO2_mumol_l", "SiOH4_mumol_l", "PO4_mumol_l"))
phys_oce.sub.phy <- subset(phys_oce, select = c("depth", "temp_deg_c", "salinity", "flurom_arbit"))

# cases need to be complete.cases AND the duplicates need to be excluded:
columns2keep.all <- complete.cases(phys_oce.sub.all) & phys_oce$keep == 1  
columns2keep.phy <- complete.cases(phys_oce.sub.phy) & phys_oce$keep == 1 
columns2keep.2014 <- columns2keep.all & phys_oce$year == 2014
columns2keep.HG <- columns2keep.all & phys_oce$latitude >= 78.5 & phys_oce$latitude <=80 & phys_oce$longitude >= -5 & phys_oce$longitude <= 11
# borders of HG after: https://www.awi.de/en/science/biosciences/deep-sea-ecology-and-technology/observatories/lter-observatory-hausgarten.html

# reduce
phys_oce.sub.2014 <- phys_oce.sub.all[columns2keep.2014,]
phys_oce.sub.HG <- phys_oce.sub.all[columns2keep.HG,]
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
sequ.sub.2014 <- sequ[columns2keep.2014,]
sequ.sub.HG <- sequ[columns2keep.HG,]

# get station names from the phys_oce datasheet:
stat_names <- subset(phys_oce, select = c("Proben_ID_intern", "date", "depth",
                                          "year", "latitude", "longitude"))
# reduce
stat_names.all <- stat_names[columns2keep.all,]
stat_names.phy <- stat_names[columns2keep.phy,]
stat_names.2014 <- stat_names[columns2keep.2014,]
stat_names.HG <- stat_names[columns2keep.HG,]

### this part was taken from the universal header! (END) ----
# remove unwanted filters, original datasets 
rm(columns2keep.all, columns2keep.phy, columns2keep.2014, columns2keep.HG, sequ, stat_names, phys_oce)
## triggers:

# calculate correlations?
calc.correl <- FALSE
# plot worldmap?
plot.wm <- FALSE


## DM --------------------------------------------------------

# physical data oder nutrients+physical data dazu
f.n0.all <- cbind(sequ.sub.all, phys_oce.sub.all)
f.n0.phy <- cbind(sequ.sub.phy, phys_oce.sub.phy)
f.n0.2014 <- cbind(sequ.sub.2014, phys_oce.sub.2014)
f.n0.HG <- cbind(sequ.sub.HG, phys_oce.sub.HG)

# remove zeros
f.n0.all <- zCompositions::cmultRepl(f.n0.all, method="CZM", label = 0)
f.n0.phy <- zCompositions::cmultRepl(f.n0.phy, method="CZM", label = 0)
f.n0.2014 <- zCompositions::cmultRepl(f.n0.2014, method="CZM", label = 0)
f.n0.HG <- zCompositions::cmultRepl(f.n0.HG, method="CZM", label = 0)

# variance-stabilizing transformation:
f.clr.all <- CoDaSeq::codaSeq.clr(f.n0.all, samples.by.row = T)
f.clr.phy <- CoDaSeq::codaSeq.clr(f.n0.phy, samples.by.row = T)
f.clr.2014 <- CoDaSeq::codaSeq.clr(f.n0.2014, samples.by.row = T)
f.clr.HG <- CoDaSeq::codaSeq.clr(f.n0.HG, samples.by.row = T)

# DM (Thilo's method)
data.all <- similarity(as.matrix(f.clr.all))
data.phy <- similarity(as.matrix(f.clr.phy))
data.2014 <- similarity(as.matrix(f.clr.2014))
data.HG <- similarity(as.matrix(f.clr.HG))

data.all <- simil_reducer(data.all)
data.phy <- simil_reducer(data.phy)
data.2014 <- simil_reducer(data.2014)
data.HG <- simil_reducer(data.HG)

lap.all <- matrixLaplacian::matrixLaplacian(data.all, plot2D = F, plot3D = F)
lap.phy <- matrixLaplacian::matrixLaplacian(data.phy, plot2D = F, plot3D = F)
lap.2014 <- matrixLaplacian::matrixLaplacian(data.2014, plot2D = F, plot3D = F)
lap.HG <- matrixLaplacian::matrixLaplacian(data.HG, plot2D = F, plot3D = F)

lap_mat.all <- lap.all$LaplacianMatrix
lap_mat.phy <- lap.phy$LaplacianMatrix
lap_mat.2014 <- lap.2014$LaplacianMatrix
lap_mat.HG <- lap.HG$LaplacianMatrix

elm.all <- eigen(lap_mat.all)
elm.phy <- eigen(lap_mat.phy)
elm.2014 <- eigen(lap_mat.2014)
elm.HG <- eigen(lap_mat.HG)

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
low.2014 <- data.frame("minus_1" = elm.2014$vectors[,ncol(data.2014)-1], "minus_2" = elm.2014$vectors[,ncol(data.2014)-2],
                       "minus_3" = elm.2014$vectors[,ncol(data.2014)-3], "minus_4" = elm.2014$vectors[,ncol(data.2014)-4],
                       "minus_5" = elm.2014$vectors[,ncol(data.2014)-5], "minus_6" = elm.2014$vectors[,ncol(data.2014)-6],
                       "minus_7" = elm.2014$vectors[,ncol(data.2014)-7], "minus_8" = elm.2014$vectors[,ncol(data.2014)-8],
                       year = as.factor(stat_names.2014$year),
                       depth = -as.numeric(stat_names.2014$depth),
                       count = 1:ncol(data.2014))
low.HG <- data.frame("minus_1" = elm.HG$vectors[,ncol(data.HG)-1], "minus_2" = elm.HG$vectors[,ncol(data.HG)-2],
                      "minus_3" = elm.HG$vectors[,ncol(data.HG)-3], "minus_4" = elm.HG$vectors[,ncol(data.HG)-4],
                      "minus_5" = elm.HG$vectors[,ncol(data.HG)-5], "minus_6" = elm.HG$vectors[,ncol(data.HG)-6],
                      "minus_7" = elm.HG$vectors[,ncol(data.HG)-7], "minus_8" = elm.HG$vectors[,ncol(data.HG)-8],
                      year = as.factor(stat_names.HG$year),
                      depth = -as.numeric(stat_names.HG$depth),
                      count = 1:ncol(data.HG))

## correlations -------
#' leaving out the first (smallest non-zero eigenvalue) and the corresponding eigenvector,
#' each eigenvector is correlated with each column in the original dataset to find correlations
#' between the eigenvectors and the predictors, determining what feature the eigenvector picks up

if (calc.correl){
  ## correlations 2014:
  cor.coef.2014 <- data.frame("eig.vec.ind" = sort(rep((length(elm.2014$values)-1):1, 10), decreasing = T),
                             "order.ev" = rep(1:10, length(elm.2014$values)-1), "p-value" = NA,
                             "best.predictor" = NA, "cor.coeff" = NA)
  coeffic <- NA
  pval <- NA
  for(i in (length(elm.2014$values)-1):1){
    for(j in 1:dim(f.clr.2014)[2]){
      coeffic[j] <- cor.test(elm.2014$vectors[, i], f.clr.2014[,j])[[4]]
      pval[j] <- cor.test(elm.2014$vectors[, i], f.clr.2014[,j])[[3]]
    }
    logi <- order(abs(coeffic), decreasing = T)[1:10]
    cor.coef.2014$best.predictor[cor.coef.2014$eig.vec.ind == i] <- paste("Var. ", colnames(f.clr.2014)[logi])
    cor.coef.2014$cor.coeff[cor.coef.2014$eig.vec.ind == i] <- coeffic[logi]
    cor.coef.2014$p.value[cor.coef.2014$eig.vec.ind == i] <- pval[logi]
  }
  
  ## correlations HAUSGARTEN:
  cor.coef.HG <- data.frame("eig.vec.ind" = sort(rep((length(elm.HG$values)-1):1, 10), decreasing = T),
                              "order.ev" = rep(1:10, length(elm.HG$values)-1), "p-value" = NA,
                              "best.predictor" = NA, "cor.coeff" = NA)
  coeffic <- NA
  pval <- NA
  for(i in (length(elm.HG$values)-1):1){
    for(j in 1:dim(f.clr.HG)[2]){
      coeffic[j] <- cor.test(elm.HG$vectors[, i], f.clr.HG[,j])[[4]]
      pval[j] <- cor.test(elm.HG$vectors[, i], f.clr.HG[,j])[[3]]
    }
    logi <- order(abs(coeffic), decreasing = T)[1:10]
    cor.coef.HG$best.predictor[cor.coef.HG$eig.vec.ind == i] <- paste("Var. ", colnames(f.clr.HG)[logi])
    cor.coef.HG$cor.coeff[cor.coef.HG$eig.vec.ind == i] <- coeffic[logi]
    cor.coef.HG$p.value[cor.coef.HG$eig.vec.ind == i] <- pval[logi]
  }
  
  
  cor.coef.HG[cor.coef.HG$order.ev == 1,]
}
#### plots ---------

ggplot(low.all, aes(y = minus_1 , x = phys_oce.sub.all$NO3_mumol_l+phys_oce.sub.all$NO2_mumol_l))+
  geom_point()+
  labs(title = "physical+nutrient subset")

ggplot(low.phy, aes(x = minus_1 , y = minus_2, col = depth))+
  geom_point()+
  labs(title = "physical subset",
       subtitle = "seems to pick up time signal!")

ggplot(low.all, aes(x = stat_names.all$longitude , y = minus_1, col = minus_1))+
  geom_point()+
  geom_smooth(method="lm", se=F)

ggplot(low.all, aes(y = stat_names.all$latitude , x = minus_1, col = minus_1))+
  geom_point()+
  geom_smooth(method="lm", se=F)


ggplot(low.all, aes(x = depth, y = minus_2, col = depth))+
  geom_point()


ggplot(low.all, aes(x = minus_2 , y = minus_3, col =depth))+
  geom_point()+
  labs(subtitle = "seems to pick up time signal!")

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