## run script for complete application of diffusion maps method incl. thresholding and analysis
rm(list=ls())

library(matrixLaplacian)
library(readxl)
## own function definitions: ---------------------------------------------------------------------
# threshapply to apply the thresholds to the data and remove all columns with no entries:
start.time <- Sys.time()

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

# standardize standardizes the data based on Thilo's recommendations (columnwise mean=0 and variance = 1)
standardize <- function(input_data){
  if(!is.matrix(input_data)) stop('input_data must be a matrix')
  
  for(i in 1:ncol(input_data)){
    input_data[, i] <- (input_data[, i] - mean(input_data[, i]))/sd(input_data[, i])
  }
  return(input_data)
}

# similarity calculates the similarity matrix, again based off Thilo's advice
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

# simil_reducer returns only the 10 most similar entries of each row
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


# set.seed(426)
# vec <- runif(440)
# testdata <- matrix(vec, ncol=11)
# 
# testdata <- standardize(testdata)
# testdata <- similarity(testdata)
# testdata <- simil_reducer(testdata)

sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";")
sequ <- t(sequ)

data <- threshapply(sequ, "90 percent")
data <- standardize(data)
data <- similarity(data)
data <- simil_reducer(data)

lap <- matrixLaplacian(data, plot2D = F, plot3D = F)
lap_mat <- lap$LaplacianMatrix

elm <- eigen(lap_mat)

write.csv2(elm$vectors, file = "DM_results_90percent.csv", row.names = F)
end.time <- Sys.time()
end.time-start.time

#----------------
## beginning with smallest eigenvalue, look at positive and negative poles of the corresponding eigenvector
# -> relate those to stations from the station data sheet

# for each column in vectors (or only for the first, say, 4) get the names of the 10 stations corresponding to
# highest and lowest values, respectively

ind.high <- matrix(NA, nrow = 10, ncol = ncol(elm$vectors))
ind.low <- matrix(NA, nrow = 10, ncol = ncol(elm$vectors))
for(i in 1:ncol(elm$vectors)){
  ind.high[,i] <- order(elm$vectors[,i], decreasing = T)[1:10] 
  ind.low[,i] <- order(elm$vectors[,i])[1:10]
}

end.time <- Sys.time()
end.time-start.time


#' namesfile[ind.high[, 1],] are the 10 most important large entries of the eigenvector corresponding to the lowest EV, 
#' namesfile[ind.low[, 1],] are the 10 most important small entries -"-
#' Thilo plotted the values of all stations in the first three EV against each other
#' that's something i will do, too, to see interesting directions of the matrix 


stat_names <- read_excel("Sequences_Hausgarten_station_data_revised.xlsx")

stat_names$year[ind.high[,3]]
stat_names$station[ind.high[,3]]

stat_names$year[ind.low[,3]]

stat_names$year[ind.high[,3]]
stat_names$year[ind.low[,2]]

# stat_names[ind.high[,3],1]
# stat_names[ind.low[,3],1]
