# script to calculate both alternative approaches to the Diffusion Map method
#' also includes the normalization part from Carina

rm(list=ls())
graphics.off()

## packages ----------------------------------------------------
library(readxl)
library(zCompositions) #' to remove non-zero entries from the raw data matrix
BiocManager::install("BiocParallel")
library(CoDaSeq) #' installed from tarball from ggloor's github repo on CoDaSeq
library(robCompositions) #' to calculate the Aitchison distance matrix
library(ggbiplot)
library(matrixLaplacian)
library(destiny)
library(rgl)
library(ggplot2)

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


## package CoDaSeq + Principal Coordinate Analysis -------------------

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";", header = F)
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")
#sequ <- data.frame(t(sequ))

stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")


# standardize data:
f.n0 <- zCompositions::cmultRepl(sequ, method="CZM", label = 0)
f.clr <- CoDaSeq::codaSeq.clr(f.n0, samples.by.row = T)


## traditional analysis ---------------------------------------------
clrdists <- robCompositions::aDist(f.clr)

#' from:
#' https://www.datacamp.com/community/tutorials/pca-analysis-r
#' to perform elegant visualisations of PCA results 

# PCA 
sequ.pca <- prcomp(clrdists, center = T, scale. = T)
ggbiplot(sequ.pca, labels = stat_names$Proben_ID_intern)

## diffusion map alt 1 ----------------------------------------------------
#' takes the f.n0 base data line and applies the diffusion maps procedure
#' to it

data <- similarity(as.matrix(f.n0))
data <- simil_reducer(data)

lap <- matrixLaplacian(data, plot2D = F, plot3D = F)
lap_mat <- lap$LaplacianMatrix

elm <- eigen(lap_mat)

ind.high <- matrix(NA, nrow = 10, ncol = ncol(elm$vectors))
ind.low <- matrix(NA, nrow = 10, ncol = ncol(elm$vectors))
for(i in 1:ncol(elm$vectors)){
  ind.high[,i] <- order(elm$vectors[,i], decreasing = T)[1:10] 
  ind.low[,i] <- order(elm$vectors[,i])[1:10]
}

## plotting results

high <- data.frame("first" = elm$vectors[,1], "second" = elm$vectors[,2],
                   year = as.factor(stat_names$year),
                   depth = -as.numeric(stat_names$depth),
                   count = 1:166)
low <- data.frame("penult" = elm$vectors[,165], "antepen" = elm$vectors[,164],
                  year = as.factor(stat_names$year),
                  depth = -as.numeric(stat_names$depth),
                  count = 1:166)

ggplot(high, aes(x = first, y = second, col = year))+geom_point()


ggplot(low, aes(x = penult, y = antepen, col = depth))+
  geom_point()+
  labs(title = "Plot of the penult against the antepen",
       subtitle = "seems to pick up depth signal!")


antepen.cor <- cor.test(low$depth, low$antepen)
antepen.lm <- lm(antepen ~ depth, data = low)

subt <- paste("correlation coefficient ", round(antepen.cor$estimate,2),
              " for antepen~depth, R^2=", round(summary(antepen.lm)$r.squared,3))
pen.cor <- cor.test(low$depth, low$penult)

  
ggplot(low, aes(y = antepen, x = depth, col = year))+
  geom_point()+
  geom_smooth(method = "lm",
              inherit.aes = FALSE,
              aes(x=depth, y = antepen),
              se=F)+
  labs(title = "depth against antepen. EV",
       subtitle = subt)



paste(stat_names$station[ind.high[,1]], stat_names$year[ind.high[,165]], stat_names$depth[ind.high[,1]], sep = "+")

## diffusion map alt 2 ------------------------------------------------------

dm <- destiny::DiffusionMap(f.n0, n_eigs = 164, k = 10)


#' old plots
# plot(dm, 1:2, col = as.factor(stat_names$year))   # colour legend correct??
# 
# plot3d(eigenvectors(dm)[, 1:3], col = as.factor(stat_names$depth), size = 3)
# plot3d(eigenvectors(dm)[, 1:3], col = unique(as.factor(stat_names$year)), size = 3)
# legend3d('topright', legend = unique(as.factor(stat_names$year)),pch = 16, col = unique(as.factor(stat_names$year)), cex=1)

dmextract <- data.frame(DC1 = dm$DC1, DC2 = dm$DC2,
                        year = as.factor(stat_names$year),
                        depth = -as.numeric(stat_names$depth))


ggplot(data = dmextract, aes(x = DC1, y = DC2, col = depth))+
  geom_point()+
  labs(title = "destiny DM",
       subtitle = "auch Tiefensignal?")

ggplot(data = dmextract, aes(x = -depth, y = DC2, col = depth))+
  geom_point()+
  geom_smooth(method = "lm",
              se=F)


destdm.cor <- cor.test(dmextract$depth, dmextract$DC2)
destdm.lm <- lm(DC2 ~ depth, data = dmextract)
summary(destdm.lm)
