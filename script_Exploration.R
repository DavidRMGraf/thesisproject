# Data Exploration Script
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
library(vegan)
library(ade4)
## script -----------------------------------------------------
stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

statdiff <- length(stat_names$Proben_ID_intern)-length(unique(stat_names$Proben_ID_intern)); statdiff
#' difference of 11 (!) stations that are duplicated in the dataset!
#' 

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)
sequ.red <- threshapply(sequ, "0.05 percent")

## script portion to count duplicates (omitted) ---------------
# counter <- 0
# for (i in 1:nrow(sequ)){
#   for (j in 1:nrow(sequ)){
#     if(i!=j){
#       if(sum(sequ[i,] == sequ[j,])==length(sequ[1,])){
#         counter <- counter+1
#       }
#     }
#   }
#   print(i)
# }
# counter

#' counter is still at 0 after this comparison,
#' so no complete duplicates exist in the dataset - which is weird, since we have 11 identical 'Proben_ID_intern' rows!
#' 

## Prinicple coordinate analysis ------------------------------------------------
#' PCoA with 'vegan' as described in GUSTA ME:
#' https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/principal-coordinates-analysis
#' 
#' PCoA is performed on two indices of dissimilarity as calculated on the full and filtered dataset, respectively

# PCoA on Morisita-Horn dissimilarity indices:
MHdist <- vegdist(sequ, method="morisita")
pcoa.1 <- cmdscale(MHdist)

MHdist2 <- vegdist(sequ.red, method="morisita")
pcoa.2 <- cmdscale(MHdist2)

# PCoA on a more classical Bray-Curtis indices:

BCdist <- vegdist(sequ, method = "bray")
pcoa.3 <- cmdscale(BCdist)

BCdist2 <- vegdist(sequ.red, method = "bray")
pcoa.4 <- cmdscale(BCdist2)

# vizualisations:
opar <- par()
par(mfrow = c(1,2))

ordiplot(pcoa.1, main = "PCoA-MH-full")
ordiplot(pcoa.2, main = "PCoA-MH-filt")

ordiplot(pcoa.3, main = "PCoA-BC-full")
ordiplot(pcoa.4, main = "PCoA-BC-filt")

par(opar)

# comparison between full and filtered datasets:
#' as described here:
#' https://sites.google.com/site/mb3gustame/hypothesis-tests/the-mantel-test

mantel(MHdist, MHdist2)
# p-value << 0.05 -> it is likely that there is a linear correlation between the two matrices
mantel(BCdist, BCdist2)
# same here

#' General Conclusion: distance matrices for two indices of distance (Bray-Curtis and Morisita-Horn) do not significantly
#' differ between the filtered and the full dataset when 0.05-percent-filtering is applied!

