# Principal Component Analysis Script
rm(list=ls())
graphics.off()
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

## packages ---------------------------------------------------

library(ggplot2)
library(readxl)
## script -----------------------------------------------------

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)

sequ <- threshapply(sequ, "0.05 percent")
colnames(sequ) <- 1:ncol(sequ)

stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")


## alternative 1 - besser    ---------

library(ggbiplot)
#' from:
#' https://www.datacamp.com/community/tutorials/pca-analysis-r
#' to perform elegant visualisations of PCA results 

HGpca.1 <- prcomp(sequ, center = T, scale. = T)

ggbiplot(HGpca.1, labels = stat_names$year)

## alternative 2 ---------------------

library(labdsv)
#' from:
#' http://ecology.msu.montana.edu/labdsv/R/labs/lab7/lab7.html
#' to perfrom an alternative approach to pca

HGpca.2 <- pca(sequ, cor = T, dim = 2)

plot(HGpca.2)
