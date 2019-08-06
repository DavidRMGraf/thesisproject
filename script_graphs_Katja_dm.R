# script to create graphs for Katja und so
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
library(gridExtra)

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
f.n0.input.all <- zCompositions::cmultRepl(f.n0.all, method="CZM", label = 0)
f.n0.input.phy <- zCompositions::cmultRepl(f.n0.phy, method="CZM", label = 0)
f.n0.input.2014 <- zCompositions::cmultRepl(f.n0.2014, method="CZM", label = 0)
f.n0.input.HG <- zCompositions::cmultRepl(f.n0.HG, method="CZM", label = 0)

# variance-stabilizing transformation:
f.clr.all <- CoDaSeq::codaSeq.clr(f.n0.input.all, samples.by.row = T)
f.clr.phy <- CoDaSeq::codaSeq.clr(f.n0.input.phy, samples.by.row = T)
f.clr.2014 <- CoDaSeq::codaSeq.clr(f.n0.input.2014, samples.by.row = T)
f.clr.HG <- CoDaSeq::codaSeq.clr(f.n0.input.HG, samples.by.row = T)

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
                      "minus_7" = elm.all$vectors[,ncol(data.all)-7], "minus_8" = elm.all$vectors[,ncol(data.all)-8])

low.phy <- data.frame("minus_1" = elm.phy$vectors[,ncol(data.phy)-1], "minus_2" = elm.phy$vectors[,ncol(data.phy)-2],
                      "minus_3" = elm.phy$vectors[,ncol(data.phy)-3], "minus_4" = elm.phy$vectors[,ncol(data.phy)-4],
                      "minus_5" = elm.phy$vectors[,ncol(data.phy)-5], "minus_6" = elm.phy$vectors[,ncol(data.phy)-6],
                      "minus_7" = elm.phy$vectors[,ncol(data.phy)-7], "minus_8" = elm.phy$vectors[,ncol(data.phy)-8])

low.2014 <- data.frame("minus_1" = elm.2014$vectors[,ncol(data.2014)-1], "minus_2" = elm.2014$vectors[,ncol(data.2014)-2],
                       "minus_3" = elm.2014$vectors[,ncol(data.2014)-3], "minus_4" = elm.2014$vectors[,ncol(data.2014)-4],
                       "minus_5" = elm.2014$vectors[,ncol(data.2014)-5], "minus_6" = elm.2014$vectors[,ncol(data.2014)-6],
                       "minus_7" = elm.2014$vectors[,ncol(data.2014)-7], "minus_8" = elm.2014$vectors[,ncol(data.2014)-8])

low.HG <- data.frame("minus_1" = elm.HG$vectors[,ncol(data.HG)-1], "minus_2" = elm.HG$vectors[,ncol(data.HG)-2],
                     "minus_3" = elm.HG$vectors[,ncol(data.HG)-3], "minus_4" = elm.HG$vectors[,ncol(data.HG)-4],
                     "minus_5" = elm.HG$vectors[,ncol(data.HG)-5], "minus_6" = elm.HG$vectors[,ncol(data.HG)-6],
                     "minus_7" = elm.HG$vectors[,ncol(data.HG)-7], "minus_8" = elm.HG$vectors[,ncol(data.HG)-8])

## graphs ---------
graphics.off

# p1 <- ggplot(data = low.2014, aes(x = minus_1, y = minus_2))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p2 <- ggplot(data = low.2014, aes(x = minus_2, y = minus_3))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p3 <- ggplot(data = low.2014, aes(x = minus_3, y = minus_4))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p4 <- ggplot(data = low.2014, aes(x = minus_4, y = minus_5))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p5 <- ggplot(data = low.2014, aes(x = minus_5, y = minus_6))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p6 <- ggplot(data = low.2014, aes(x = minus_6, y = minus_7))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p7 <- ggplot(data = low.2014, aes(x = minus_7, y = minus_8))+
#   geom_point()+
#   labs(title = "subset 2014")+
#   coord_fixed(ratio = 1)
# 
# p11 <- ggplot(data = low.all, aes(x = minus_1, y = minus_2))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p12 <- ggplot(data = low.all, aes(x = minus_2, y = minus_3))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p13 <- ggplot(data = low.all, aes(x = minus_3, y = minus_4))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p14 <- ggplot(data = low.all, aes(x = minus_4, y = minus_5))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p15 <- ggplot(data = low.all, aes(x = minus_5, y = minus_6))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p16 <- ggplot(data = low.all, aes(x = minus_6, y = minus_7))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p17 <- ggplot(data = low.all, aes(x = minus_7, y = minus_8))+
#   geom_point()+
#   labs(title = "subset all")+
#   coord_fixed(ratio = 1)
# 
# p21 <- ggplot(data = low.HG, aes(x = minus_1, y = minus_2))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p22 <- ggplot(data = low.HG, aes(x = minus_2, y = minus_3))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p23 <- ggplot(data = low.HG, aes(x = minus_3, y = minus_4))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p24 <- ggplot(data = low.HG, aes(x = minus_4, y = minus_5))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p25 <- ggplot(data = low.HG, aes(x = minus_5, y = minus_6))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p26 <- ggplot(data = low.HG, aes(x = minus_6, y = minus_7))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p27 <- ggplot(data = low.HG, aes(x = minus_7, y = minus_8))+
#   geom_point()+
#   labs(title = "subset HG")+
#   coord_fixed(ratio = 1)
# 
# p31 <- ggplot(data = low.phy, aes(x = minus_1, y = minus_2))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p32 <- ggplot(data = low.phy, aes(x = minus_2, y = minus_3))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p33 <- ggplot(data = low.phy, aes(x = minus_3, y = minus_4))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p34 <- ggplot(data = low.phy, aes(x = minus_4, y = minus_5))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p35 <- ggplot(data = low.phy, aes(x = minus_5, y = minus_6))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p36 <- ggplot(data = low.phy, aes(x = minus_6, y = minus_7))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)
# 
# p37 <- ggplot(data = low.phy, aes(x = minus_7, y = minus_8))+
#   geom_point()+
#   labs(title = "subset phy")+
#   coord_fixed(ratio = 1)

p1 <- ggplot(data = low.2014, aes(x = minus_1, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")

p2 <- ggplot(data = low.2014, aes(x = minus_2, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")

p3 <- ggplot(data = low.2014, aes(x = minus_3, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p4 <- ggplot(data = low.2014, aes(x = minus_4, y = minus_5))+
  geom_point()+
  labs(title = "subset 2014")

p5 <- ggplot(data = low.2014, aes(x = minus_5, y = minus_6))+
  geom_point()+
  labs(title = "subset 2014")

p6 <- ggplot(data = low.2014, aes(x = minus_6, y = minus_7))+
  geom_point()+
  labs(title = "subset 2014")

p7 <- ggplot(data = low.2014, aes(x = minus_7, y = minus_8))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = low.all, aes(x = minus_1, y = minus_2))+
  geom_point()+
  labs(title = "subset all")

p12 <- ggplot(data = low.all, aes(x = minus_2, y = minus_3))+
  geom_point()+
  labs(title = "subset all")

p13 <- ggplot(data = low.all, aes(x = minus_3, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p14 <- ggplot(data = low.all, aes(x = minus_4, y = minus_5))+
  geom_point()+
  labs(title = "subset all")

p15 <- ggplot(data = low.all, aes(x = minus_5, y = minus_6))+
  geom_point()+
  labs(title = "subset all")

p16 <- ggplot(data = low.all, aes(x = minus_6, y = minus_7))+
  geom_point()+
  labs(title = "subset all")

p17 <- ggplot(data = low.all, aes(x = minus_7, y = minus_8))+
  geom_point()+
  labs(title = "subset all")

p21 <- ggplot(data = low.HG, aes(x = minus_1, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")

p22 <- ggplot(data = low.HG, aes(x = minus_2, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")

p23 <- ggplot(data = low.HG, aes(x = minus_3, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

p24 <- ggplot(data = low.HG, aes(x = minus_4, y = minus_5))+
  geom_point()+
  labs(title = "subset HG")

p25 <- ggplot(data = low.HG, aes(x = minus_5, y = minus_6))+
  geom_point()+
  labs(title = "subset HG")

p26 <- ggplot(data = low.HG, aes(x = minus_6, y = minus_7))+
  geom_point()+
  labs(title = "subset HG")

p27 <- ggplot(data = low.HG, aes(x = minus_7, y = minus_8))+
  geom_point()+
  labs(title = "subset HG")

p31 <- ggplot(data = low.phy, aes(x = minus_1, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")

p32 <- ggplot(data = low.phy, aes(x = minus_2, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")

p33 <- ggplot(data = low.phy, aes(x = minus_3, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p34 <- ggplot(data = low.phy, aes(x = minus_4, y = minus_5))+
  geom_point()+
  labs(title = "subset phy")

p35 <- ggplot(data = low.phy, aes(x = minus_5, y = minus_6))+
  geom_point()+
  labs(title = "subset phy")

p36 <- ggplot(data = low.phy, aes(x = minus_6, y = minus_7))+
  geom_point()+
  labs(title = "subset phy")

p37 <- ggplot(data = low.phy, aes(x = minus_7, y = minus_8))+
  geom_point()+
  labs(title = "subset phy")

grid.arrange(p11, p12, p13, p14, p15, p16, p17,
             p1, p2, p3, p4, p5, p6, p7,
             p21, p22, p23, p24, p25, p26, p27,
             p31, p32, p33, p34, p35, p36, p37, nrow = 4)
graphics.off()

phys_oce.sub.2014$nitrogen_species <- phys_oce.sub.2014$NO3_mumol_l+phys_oce.sub.2014$NO2_mumol_l
phys_oce.sub.HG$nitrogen_species <- phys_oce.sub.HG$NO3_mumol_l+phys_oce.sub.HG$NO2_mumol_l
phys_oce.sub.all$nitrogen_species <- phys_oce.sub.all$NO3_mumol_l+phys_oce.sub.all$NO2_mumol_l

tot_sub.2014 <- cbind(low.2014, phys_oce.sub.2014)
tot_sub.HG <- cbind(low.HG, phys_oce.sub.HG)
tot_sub.all <- cbind(low.all, phys_oce.sub.all)
tot_sub.phy <- cbind(low.phy, phys_oce.sub.phy)

p1 <- ggplot(data = tot_sub.2014, aes(x = depth, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = depth, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = depth, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = depth, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = depth, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = depth, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = depth, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = depth, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = depth, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = depth, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = depth, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = depth, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = depth, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = depth, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = depth, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = depth, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = temp_deg_c, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = temp_deg_c, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = temp_deg_c, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = temp_deg_c, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = temp_deg_c, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = temp_deg_c, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = temp_deg_c, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = temp_deg_c, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = temp_deg_c, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = temp_deg_c, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = temp_deg_c, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = temp_deg_c, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = temp_deg_c, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = temp_deg_c, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = temp_deg_c, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = temp_deg_c, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = salinity, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = salinity, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = salinity, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = salinity, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = salinity, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = salinity, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = salinity, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = salinity, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = salinity, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = salinity, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = salinity, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = salinity, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = salinity, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = salinity, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = salinity, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = salinity, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = flurom_arbit, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = flurom_arbit, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = flurom_arbit, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = flurom_arbit, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = flurom_arbit, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = flurom_arbit, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = flurom_arbit, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = flurom_arbit, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = flurom_arbit, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = flurom_arbit, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = flurom_arbit, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = flurom_arbit, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = flurom_arbit, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = flurom_arbit, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = flurom_arbit, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = flurom_arbit, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = SiOH4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = SiOH4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = SiOH4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = SiOH4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p21 <- ggplot(data = tot_sub.all, aes(x = SiOH4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = SiOH4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = SiOH4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = SiOH4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = SiOH4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = SiOH4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = SiOH4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = SiOH4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34, nrow = 3)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = PO4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = PO4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = PO4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = PO4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p21 <- ggplot(data = tot_sub.all, aes(x = PO4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = PO4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = PO4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = PO4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = PO4_mumol_l, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = PO4_mumol_l, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = PO4_mumol_l, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = PO4_mumol_l, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34, nrow = 3)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = nitrogen_species, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = nitrogen_species, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = nitrogen_species, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = nitrogen_species, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p21 <- ggplot(data = tot_sub.all, aes(x = nitrogen_species, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = nitrogen_species, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = nitrogen_species, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = nitrogen_species, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = nitrogen_species, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = nitrogen_species, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = nitrogen_species, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = nitrogen_species, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34, nrow = 3)

graphics.off()

tot_sub.2014 <- cbind(tot_sub.2014, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
tot_sub.all <- cbind(tot_sub.all, latitude = stat_names.all$latitude, longitude = stat_names.all$longitude)
tot_sub.phy <- cbind(tot_sub.phy, latitude = stat_names.phy$latitude, longitude = stat_names.phy$longitude)
tot_sub.HG <- cbind(tot_sub.HG, latitude = stat_names.HG$latitude, longitude = stat_names.HG$longitude)


p1 <- ggplot(data = tot_sub.2014, aes(x = latitude, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = latitude, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = latitude, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = latitude, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = latitude, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = latitude, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = latitude, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = latitude, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = latitude, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = latitude, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = latitude, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = latitude, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = latitude, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = latitude, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = latitude, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = latitude, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()

p1 <- ggplot(data = tot_sub.2014, aes(x = longitude, y = minus_1))+
  geom_point()+
  labs(title = "subset 2014")
p2 <- ggplot(data = tot_sub.2014, aes(x = longitude, y = minus_2))+
  geom_point()+
  labs(title = "subset 2014")
p3 <- ggplot(data = tot_sub.2014, aes(x = longitude, y = minus_3))+
  geom_point()+
  labs(title = "subset 2014")
p4 <- ggplot(data = tot_sub.2014, aes(x = longitude, y = minus_4))+
  geom_point()+
  labs(title = "subset 2014")

p11 <- ggplot(data = tot_sub.phy, aes(x = longitude, y = minus_1))+
  geom_point()+
  labs(title = "subset phy")
p12 <- ggplot(data = tot_sub.phy, aes(x = longitude, y = minus_2))+
  geom_point()+
  labs(title = "subset phy")
p13 <- ggplot(data = tot_sub.phy, aes(x = longitude, y = minus_3))+
  geom_point()+
  labs(title = "subset phy")
p14 <- ggplot(data = tot_sub.phy, aes(x = longitude, y = minus_4))+
  geom_point()+
  labs(title = "subset phy")

p21 <- ggplot(data = tot_sub.all, aes(x = longitude, y = minus_1))+
  geom_point()+
  labs(title = "subset all")
p22 <- ggplot(data = tot_sub.all, aes(x = longitude, y = minus_2))+
  geom_point()+
  labs(title = "subset all")
p23 <- ggplot(data = tot_sub.all, aes(x = longitude, y = minus_3))+
  geom_point()+
  labs(title = "subset all")
p24 <- ggplot(data = tot_sub.all, aes(x = longitude, y = minus_4))+
  geom_point()+
  labs(title = "subset all")

p31 <- ggplot(data = tot_sub.HG, aes(x = longitude, y = minus_1))+
  geom_point()+
  labs(title = "subset HG")
p32 <- ggplot(data = tot_sub.HG, aes(x = longitude, y = minus_2))+
  geom_point()+
  labs(title = "subset HG")
p33 <- ggplot(data = tot_sub.HG, aes(x = longitude, y = minus_3))+
  geom_point()+
  labs(title = "subset HG")
p34 <- ggplot(data = tot_sub.HG, aes(x = longitude, y = minus_4))+
  geom_point()+
  labs(title = "subset HG")

grid.arrange(p21, p22, p23, p24,
             p1, p2, p3, p4,
             p31, p32, p33, p34,
             p11, p12, p13, p14, nrow = 4)

graphics.off()