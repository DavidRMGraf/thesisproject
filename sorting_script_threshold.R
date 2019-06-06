rm(list=ls())

set.seed(426)
vec <- runif(440)
mat <- matrix(vec, ncol=11)
rm(vec)

## 90%-threshold
for (i in 1:ncol(mat)){
  mat[order(mat[, i], decreasing=T)[cumsum(mat[order(mat[, i], decreasing=T), i])/sum(mat[, i])>0.90], i] <- 0
}
# 
# 
# ## 95%-threshold
# for (i in 1:ncol(mat)){
#   mat[order(mat[, i], decreasing=T)[cumsum(mat[order(mat[, i], decreasing=T), i])/sum(mat[, i])>0.95], i] <- 0
# }
# 
# ## 99%-threshold
# for (i in 1:ncol(mat)){
#   mat[order(mat[, i], decreasing=T)[cumsum(mat[order(mat[, i], decreasing=T), i])/sum(mat[, i])>0.99], i] <- 0
# }
# 
