# variance normalization application and Principa Component Analysis script
rm(list=ls())
graphics.off()

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


## package CoDaSeq + Principal Coordinate Analysis -------------------

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")
#sequ <- data.frame(t(sequ))

stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

library(zCompositions) #' to remove non-zero entries from the raw data matrix
library(CoDaSeq) #' installed from tarball from ggloor's github repo on CoDaSeq
library(robCompositions) #' to calculate the Aitchison distance matrix
library(ggbiplot)
library(matrixLaplacian)

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



paste(stat_names$station[ind.high[,1]], stat_names$year[ind.high[,165]], stat_names$depth[ind.high[,1]], sep = "+")
# paste(stat_names$station[ind.high[,2]], stat_names$year[ind.high[,2]], stat_names$depth[ind.high[,2]], sep = "+")
# paste(stat_names$station[ind.high[,3]], stat_names$year[ind.high[,3]], stat_names$depth[ind.high[,3]], sep = "+")

## diffusion map alt 2 ------------------------------------------------------
library(destiny)
library(rgl)
library(ggplot2)

dm <- destiny::DiffusionMap(f.n0, n_eigs = 164)
plot(dm, 1:2, col = as.factor(stat_names$year))   # colour legend correct??

plot3d(eigenvectors(dm)[, 1:3], col = as.factor(stat_names$depth), size = 3)
plot3d(eigenvectors(dm)[, 1:3], col = unique(as.factor(stat_names$year)), size = 3)
legend3d('topright', legend = unique(as.factor(stat_names$year)),pch = 16, col = unique(as.factor(stat_names$year)), cex=1)

qplot(DC1, DC2, data = dm, colour =as.factor(stat_names$depth)) +scale_color_cube_helix()
qplot(DC1, DC2, data = dm, colour =as.factor(stat_names$year)) +scale_color_cube_helix()
qplot(DC1, DC2, data = dm, colour =as.factor(stat_names$latitude)) +scale_color_cube_helix()
qplot(DC1, DC2, data = dm, colour =as.factor(stat_names$longitude)) +scale_color_cube_helix()

qplot(y = eigenvalues(dm))+
  theme_minimal()+
  labs(x ='Diffusion component (DC)', y ='Eigenvalue', title = "dm")

qplot(y = eigen(lap_mat, only.values = T)$values)+
  theme_minimal()+
  labs(x ='Diffusion component (DC)', y ='Eigenvalue', title = "zu Fuß")

# ## package DESeq [obsolete] -------------------
# 
# # biocLite("DESeq")
# # library(DESeq2)
# # 
# # sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
# # sequ <- t(sequ)
# # sequ <- threshapply(sequ, "0.05 percent")
# # #sequ <- data.frame(t(sequ))
# # 
# # stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")
# # 
# # 
# # # sequ <- newCountDataSet(sequ, conditions = stat_names$station)
# # # cell_ij in the i-th row and the j-th column includes reads of gene i in sample j
# # # this is how it is supposed to be, based on the package documentation (https://fckaf.de/zqn)
# # 
# # sequ <- estimateSizeFactors(sequ)
# # 
# # sequBlind <- estimateDispersions(sequ, method = "blind")
# # 
# # # der resultiert mit jeder Möglichen Method in einem Fehler weil
# # # 'glm.fit: algorithm did not converge'
# # 
# # vsd = varianceStabilizingTransformation(sequ)
# # 
# # ?estimateDispersions
# # 
# 
