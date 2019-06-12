rm(list=ls())

# any random matrix for functionality testing:
set.seed(426)
vec <- runif(44000)
rand.mat <- matrix(vec, ncol=110)
rm(vec)

# system time at start:
#start.time <- Sys.time()

### description ----------------------------
# virtual dataframe is structured with each species in it's own column with location as rows

### per-station species thresholds: ----------------
## only keep entries if they are part of the lowest number of entries for each OTU that constitute 90/95/99 percent of their OTU

# threshapply <- function(input_data, method){
#   if(!is.character(method)) stop('method must be character')
#   if(!is.matrix(input_data)) stop('input_data must be a matrix')
#   if (method=="90 percent"){
#     for (i in 1:nrow(input_data)){
#       output_data <- input_data
#       output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.90]] <- 0
#     }
#   }else if(method == "95 percent"){
#     for (i in 1:nrow(input_data)){
#       output_data <- input_data
#       output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.95]] <- 0
#     }
#   }else if(method == "99 percent"){
#     for (i in 1:nrow(input_data)){
#       output_data <- input_data
#       output_data[i, order(input_data[i, ], decreasing=T)[cumsum(input_data[i, order(input_data[i, ], decreasing=T)])/sum(input_data[i, ])>0.99]] <- 0
#     }
#   }else if(method == "0.05 percent"){
#     keep.cols <- colSums(input_data)/sum(input_data)>=5e-04
#     output_data <- input_data[, keep.cols]
#   }else if(method == "0.005 percent"){
#     keep.cols <- colSums(input_data)/sum(input_data)>=5e-05
#     output_data <- input_data[, keep.cols]
#   }
#   output_data <- output_data[, colSums(output_data)!=0]
#   return(output_data)
# }


## 90%-threshold 
# for (i in 1:nrow(rand.mat)){
#   rand.mat[i, order(rand.mat[i, ], decreasing=T)[cumsum(rand.mat[i, order(rand.mat[i, ], decreasing=T)])/sum(rand.mat[i, ])>0.90]] <- 0
# }

## 95%-threshold
# for (i in 1:nrow(rand.mat)){
#   rand.mat[i, order(rand.mat[i, ], decreasing=T)[cumsum(rand.mat[i, order(rand.mat[i, ], decreasing=T)])/sum(rand.mat[i, ])>0.95]] <- 0
# }

## 99%-threshold
# for (i in 1:nrow(rand.mat)){
#   rand.mat[i, order(rand.mat[i, ], decreasing=T)[cumsum(rand.mat[i, order(rand.mat[i, ], decreasing=T)])/sum(rand.mat[i, ])>0.99]] <- 0
# }

### total overall threshold: ---------------
## only keep OTUs that contribute at least 0.05% of all individuals in the analysis 
# keep.cols <- colSums(rand.mat)/sum(rand.mat)>=5e-04
# rand.mat <- rand.mat[, keep.cols]

## or 0.005%, respectively:
# keep.cols <- colSums(rand.mat)/sum(rand.mat)>=5e-05
# rand.mat <- rand.mat[, keep.cols]

### sparse Matrix format to save space -----
## to create 'sparse' matrices

# library(Matrix) 

## http://www.johnmyleswhite.com/notebook/2011/10/31/using-sparse-matrices-in-r/

### system time at end and difference between start and end (RUNTIME) ------------
# end.time <- Sys.time()
# end.time - start.time

## old (wrong) threshold code ----------------------------------------------------
# for (i in 1:ncol(rand.mat)){
#   rand.mat[order(rand.mat[, i], decreasing=T)[cumsum(rand.mat[order(rand.mat[, i], decreasing=T), i])/sum(rand.mat[, i])>0.90], i] <- 0
# } 