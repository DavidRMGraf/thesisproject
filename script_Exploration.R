# Data Exploration Script
rm(list=ls())
graphics.off()

## packages ---------------------------------------------------
library(ggplot2)
library(readxl)
library(vegan)
## script -----------------------------------------------------
stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

statdiff <- length(stat_names$Proben_ID_intern)-length(unique(stat_names$Proben_ID_intern)); statdiff
#' difference of 11 (!) stations that are duplicated in the dataset!
#' 

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)

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

## PCoA ------------------------------------------------
#' Prinicple coordinate analysis with 'vegan' as described in
#' Gusta ME's https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/principal-coordinates-analysis
#' PCoA performed on a morisita-horn dissimilarity matrix:

MHdist <- vegdist(sequ, method="morisita")
pcoa.1 <- cmdscale(MHdist)
ordiplot(pcoa.1)

#' PCoA on a more classical Bray-Curtis matrix:

BCdist <- vegdist(sequ, method = "bray")
pcoa.2 <- cmdscale(BCdist)
ordiplot(pcoa.2)



