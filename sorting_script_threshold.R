rm(list=ls())

set.seed(426)
vec <- runif(440)
mat <- matrix(vec, ncol=11)
rm(vec)

for (i in 1:ncol(mat)){
  mat[order(mat[, i], decreasing=T)[cumsum(mat[order(mat[, i], decreasing=T), i])/sum(mat[, i])>0.95], i] <- 0
}
