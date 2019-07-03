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

f.n0 <- zCompositions::cmultRepl(sequ, method="CZM", label = 0)
f.clr <- CoDaSeq::codaSeq.clr(f.n0, samples.by.row = T)
clrdists <- robCompositions::aDist(f.clr)

#' from:
#' https://www.datacamp.com/community/tutorials/pca-analysis-r
#' to perform elegant visualisations of PCA results 

# PCA ---------------
sequ.pca <- prcomp(clrdists, center = T, scale. = T)
ggbiplot(sequ.pca, labels = stat_names$Proben_ID_intern)

# ## package DESeq -------------------
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
