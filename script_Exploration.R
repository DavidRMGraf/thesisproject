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

# PCoA based on Chao-matrix:

Chdist <- vegdist(sequ, method = "chao")
ord.5 <- cmdscale(Chdist)

Chdist2 <- vegdist(sequ.red, method = "chao")
ord.6 <- cmdscale(Chdist2)

# PCoA on Morisita-Horn dissimilarity indices:
# MHdist <- vegdist(sequ, method="morisita")
# ord.1 <- cmdscale(MHdist)
# 
# MHdist2 <- vegdist(sequ.red, method="morisita")
# ord.2 <- cmdscale(MHdist2)
# 
# # PCoA on a more classical Bray-Curtis indices:
# 
# BCdist <- vegdist(sequ, method = "bray")
# ord.3 <- cmdscale(BCdist)
# 
# BCdist2 <- vegdist(sequ.red, method = "bray")
# ord.4 <- cmdscale(BCdist2)
#
# # PCoA based on Euclidean-matrix:
# 
# Eudist <- vegdist(sequ, method = "euclidean")
# ord.7 <- cmdscale(Eudist)
# 
# Eudist2 <- vegdist(sequ.red, method = "euclidean")
# ord.8 <- cmdscale(Eudist2)
# 
# # PCoA based on Jaccard index
#
# Jadist <- vegdist(sequ, method = "jaccard")
# ord.9 <- cmdscale(Jadist)
# 
# Jadist2 <- vegdist(sequ.red, method = "jaccard")
# ord.10 <- cmdscale(Jadist2)

# vizualisations:
opar <- par()
par(mfrow = c(1,2))

#ordiplot(ord.5, main = "PCoA-Ch-full")
#ordiplot(ord.6, main = "PCoA-Ch-filt")
#
# ordiplot(ord.1, main = "PCoA-MH-full")
# ordiplot(ord.2, main = "PCoA-MH-filt")
# 
# ordiplot(ord.3, main = "PCoA-BC-full")
# ordiplot(ord.4, main = "PCoA-BC-filt")
# 
# ordiplot(ord.7, main = "PCoa-Eu-full")
# ordiplot(ord.8, main = "PCoa-Eu-filt")
# 
# ordiplot(ord.9, main = "PCoa-Ja-full")
# ordiplot(ord.10, main = "PCoa-Ja-filt")
#' -> clear difference between filtered and non-filtered dataset as Jaccard index
#' is based on presence/absence data


par(opar)

# comparison between full and filtered datasets:
#' as described here:
#' https://sites.google.com/site/mb3gustame/hypothesis-tests/the-mantel-test

#mantel(MHdist, MHdist2)
# p-value << 0.05 -> it is likely that there is a linear correlation between the two matrices
#mantel(BCdist, BCdist2)
# same here

mantel(Chdist, Chdist2)
# p-value < 0.05 -> linear correlation between the two distance matrices

# mantel(Eudist, Eudist2)

#' General Conclusion: distance matrices for two indices of distance (Bray-Curtis and Morisita-Horn) do not significantly
#' differ between the filtered and the full dataset when 0.05-percent-filtering is applied!

## clustering -----------------------------------------
clust.Ch <- hclust(Chdist2, method = "complete")
cutoff <- cutree(clust.Ch, 3)

ordiplot(Chdist2)
ordihull(Chdist2, cutoff)
ordispider(Chdist2, cutoff, col = 'lightblue')
ordiellipse(Chdist2, cutoff)

#' wenn ich mir das hier so ansehe scheint clustering bei den OTU-Abundanzen recht sinnlos!

## data transformation nach Carin ---------------------
#' codaSeq-package is not available for R version 3.3.3!
#' alternatives?

# library(zCompositions)
# library(codaSeq)
# 
# f.n0 <- zCompositions::cmultRepl(sequ, method="CZM", label=0)
# f.clr <- codaSeq::codaSeq.clr(f.n0)