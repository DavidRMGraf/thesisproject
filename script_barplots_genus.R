## barplot viz script of genus/taxon/etc level abundances
rm(list=ls())
graphics.off()
## packages --------------
library(tidyverse)
library(plyr)
## functions -------------
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
## script ------------------
#' data plumbing:

otuabu <- read.csv("Sequences_Hausgarten2009-2016+names.csv", sep = ";", header = F)
otuabu <- otuabu[-c(1:5), ]
otuabu[,1] <- as.numeric(as.vector.factor(otuabu[,1]))
rn <- otuabu[, 1]
rm(list = c("otuabu"))

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";", header = F)
sequ <- t(sequ)
colnames(sequ) <- rn

#' apply threshold to data - remember to adapt the threshold used for here after threshold in the DM/PCA is chosen!
sequ.red <- threshapply(sequ, "0.05 percent")
abun_otu <- colSums(sequ.red)
otu_names <- colnames(sequ.red)
raw_bar_dat <- data.frame(count = colSums(sequ.red), id = as.character(colnames(sequ.red)), stringsAsFactors = F)

taxo <- read.csv("Taxonomy.csv", sep = ";", header = F)
taxo$id <- as.character(rn)
colnames(taxo) <- c("domain", "clade_line", "kingdom", "phylum", "subphylum_order", "class_family", "genus", "species", "id")

taxos_left <- dplyr::inner_join(taxo, raw_bar_dat)

countdat <- ddply(taxos_left, ("subphylum_order"), summarize,
                     tot_sum = sum(count, na.rm =F))
countdat$subphylum_order <- as.character(countdat$subphylum_order)

ggplot(countdat, aes(x=subphylum_order, y = tot_sum))+
  geom_histogram(stat = "identity")+
  coord_flip()

countdat <- ddply(taxos_left, ("phylum"), summarize,
                  tot_sum = sum(count, na.rm =F))
countdat$phylum <- as.character(countdat$phylum)

ggplot(countdat, aes(x=phylum, y = tot_sum))+
  geom_histogram(stat = "identity")+
  coord_flip()
