# comparison of thresholds - script
rm(list=ls())
## reading in data --------------------------------
sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)
# data frame is now structured with OTUs in columsn and stations in rows

## thresholding variant 1----------------------------

# 90%:
sequ90 <- sequ
for (i in 1:nrow(sequ)){
  sequ90[i, order(sequ[i, ], decreasing=T)[cumsum(sequ[i, order(sequ[i, ], decreasing=T)])/sum(sequ[i, ])>0.90]] <- 0
}

sequ90 <- sequ90[, colSums(sequ90)!=0]

# 95%:
sequ95 <- sequ
for (i in 1:nrow(sequ)){
  sequ95[i, order(sequ[i, ], decreasing=T)[cumsum(sequ[i, order(sequ[i, ], decreasing=T)])/sum(sequ[i, ])>0.95]] <- 0
}

sequ95 <- sequ95[, colSums(sequ95)!=0]

# 99%:
sequ99 <- sequ
for (i in 1:nrow(sequ)){
  sequ99[i, order(sequ[i, ], decreasing=T)[cumsum(sequ[i, order(sequ[i, ], decreasing=T)])/sum(sequ[i, ])>0.99]] <- 0
}

sequ99 <- sequ99[, colSums(sequ99)!=0]

## thresholding variant 2 --------------------------
keep.cols <- colSums(sequ)/sum(sequ)>=5e-04
sequK1 <- sequ[, keep.cols]
sequK1 <- sequK1[, colSums(sequK1)!=0]

keep.cols <- colSums(sequ)/sum(sequ)>=5e-05
sequK2 <- sequ[, keep.cols]
sequK2 <- sequK2[, colSums(sequK2)!=0]

## comparison --------------------------------------
print(ncol(sequ90)) # 1255
print(ncol(sequ95)) # 2208
print(ncol(sequ99)) # 10983
print(ncol(sequK1)) # 225
print(ncol(sequK2)) # 803

## rank-abundance plots ----------------------------
# for (i in 1:nrow(sequK2)){
#   plot(1:ncol(sequK2), sort(sequK2[i,], decreasing=T))  
# }
plot(1:ncol(sequK2), sort(sequK2[7,], decreasing=T), type = "l")
plot(1:ncol(sequ90), sort(sequ90[7,], decreasing=T))
