rm(list=ls())
# to create 'sparse' matrices
library(Matrix) 
# http://www.johnmyleswhite.com/notebook/2011/10/31/using-sparse-matrices-in-r/

library(readr)

sequ <- as.data.frame(read_delim("Sequences_Hausgarten2009-2016+names.csv", 
                   ";",
                   escape_double = FALSE,
                   trim_ws = TRUE
                   ))
sequ <- sequ[-(1:4), -1]

results <- vector(length = ncol(sequ))

for(i in 1:ncol(sequ)){
  results[i] <- length(as.numeric(sequ[,i])[as.numeric(sequ[,i])>=mean(as.numeric(sequ[,i]))+3*sd(as.numeric(sequ[,i]))])
}

hist(results, main = sum(results))
