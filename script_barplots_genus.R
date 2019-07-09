## barplot viz script
rm(list=ls())
  
otuabu <- read.csv("Sequences_Hausgarten2009-2016+names.csv", sep = ";", header = F)
otuabu <- otuabu[-c(1:5), ]
otuabu[,1] <- as.numeric(as.vector.factor(otuabu[,1]))
rn <- otuabu[, 1]
rm(list = c("otuabu"))

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";", header = F)
sequ <- t(sequ)

