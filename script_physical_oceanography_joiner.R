## script to extract data for the analyses:
#' one list with all the data, excluding 2016, in a list, with levels
#' sequences, names, physical data
#' 
#' second list with ALL the data, i.e. sequences including 2016
#' but without nutrients (unavailable as of now!)

rm(list = ls())
graphics.off()

## savenames ---------------
phys.oce.export.name <- "physical_oceanography_data_all_years.rds"
sequence.export.name <- "sequences_thresh_applied.rds"

## packages ----------------
library(readxl)

## functions ---------------
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
sequ <- read.csv("Sequences_Hausgarten2009-2016_ohne_header.csv", sep = ";", header = F)
sequ <- t(sequ)
sequ <- threshapply(sequ, "0.05 percent")

saveRDS(sequ, sequence.export.name)

# reading in the data on physical oceanography
stat_names <- read.csv("stat_names_physical_data.csv", header = T)
stat_names$Proben_ID_intern <- as.character.factor(stat_names$Proben_ID_intern)
stat_names$sample_ID <- as.character.factor(stat_names$sample_ID)
stat_names$date <- as.character.factor(stat_names$date)
stat_names$ARK_cruise_ID <- as.character.factor(stat_names$ARK_cruise_ID)
stat_names$cruise <- as.character.factor(stat_names$cruise)
stat_names$Date_char <- as.character.factor(stat_names$Date_char)
stat_names$Proben_ID_intern[duplicated(stat_names$Proben_ID_intern)] <- paste0(stat_names$Proben_ID_intern[duplicated(stat_names$Proben_ID_intern)], "_secmeas")

# preparing the data on physical oceanography
phys_oce <- stat_names[, c(1,3,4,5,9,10)]
# remove all non-complete cases, leaving out the coordinates
for(i in 12:ncol(stat_names)-1){
  if(complete.cases(t(stat_names[,i]))){
    n <- names(phys_oce)
    phys_oce <- cbind(phys_oce, stat_names[,i])
    names(phys_oce) <- c(n, as.character(names(stat_names)[i]))
  }
}

# reading in the data on nutrients
stat_names <- read.csv("stat_names_nutrient_data.csv", header = T)

# joining the datasets on physical oceanography and nutrient data:
n <- names(phys_oce)
phys_oce <- cbind(phys_oce, stat_names[,13:16])
str(phys_oce)

phys_oce$Date_char <- as.character.factor(phys_oce$Date_char)

saveRDS(phys_oce, phys.oce.export.name)
