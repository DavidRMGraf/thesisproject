rm(list=ls())

# any random matrix for functionality testing:
set.seed(426)
vec <- runif(440)
rand.mat <- matrix(vec, ncol=11)
rm(vec)

# system time at start:
start.time <- Sys.time()

# description ----------------------------
# virtual dataframe is structured with each species in it's own column with location as rows

# per-species thresholds: ----------------
## only keep entries if they are part of the lowest number of entries for each OTU that constitute 90/95/99 percent of their OTU

## 90%-threshold 
for (i in 1:ncol(rand.mat)){
  rand.mat[order(rand.mat[, i], decreasing=T)[cumsum(rand.mat[order(rand.mat[, i], decreasing=T), i])/sum(rand.mat[, i])>0.90], i] <- 0
} 

## 95%-threshold
# for (i in 1:ncol(rand.mat)){
#   rand.mat[order(rand.mat[, i], decreasing=T)[cumsum(rand.mat[order(rand.mat[, i], decreasing=T), i])/sum(rand.mat[, i])>0.95], i] <- 0
# }

## 99%-threshold
# for (i in 1:ncol(rand.mat)){
#   rand.mat[order(rand.mat[, i], decreasing=T)[cumsum(rand.mat[order(rand.mat[, i], decreasing=T), i])/sum(rand.mat[, i])>0.99], i] <- 0
# }
# 

# total overall threshold: ---------------
## only keep OTUs that contribute at least 0.05% of all individuals in the analysis 

# rand.mat <- rand.mat[, -colSums(rand.mat)/sum(rand.mat)<5e-04]

## or 0.005%, respectively:

# rand.mat <- rand.mat[, -colSums(rand.mat)/sum(rand.mat)<5e-05]

# system time at end and difference between start and end (RUNTIME) ------------
end.time <- Sys.time()
end.time - start.time