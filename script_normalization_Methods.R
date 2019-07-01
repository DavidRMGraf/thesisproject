# variance normalization application script
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

## package edgeR -------------------

#' first installation stuff:
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")

#' library(edgeR)
#' 
#' ##
#' sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
#' #' no need to transform the data as the edgeR package needs rows as genes and columns as libraries!
#' 
#' # transform count data sheet to DGElist-file
#' sequ <- DGEList(counts = sequ)
#' 

## package DESeq -------------------

# biocLite("DESeq")
library(DESeq)

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")
sequ <- data.frame(t(sequ))

stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

sequ <- newCountDataSet(sequ, conditions = stat_names$station)
# cell_ij in the i-th row and the j-th column includes reads of gene i in sample j
# this is how it is supposed to be, based on the package documentation (https://fckaf.de/zqn)

sequ <- estimateSizeFactors(sequ)

sequBlind <- estimateDispersions(sequ, method = "pooled")

# der resultiert mit jeder Möglichen Method in einem Fehler weil
# 'glm.fit: algorithm did not converge'

vsd = varianceStabilizingTransformation(sequ)

?estimateDispersions


