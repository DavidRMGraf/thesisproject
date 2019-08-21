## corrplots script
rm(list=ls())
graphics.off()

setwd("~/Studium/19SS/BA/thesisproject")
# get analysis data:
source("script_analysis_all_data.r")

rm(list = ls()[-(c(29:31, 38:40))])


corrdat.2014 <- cor(cbind(Depth = phys_oce.sub.2014$depth,
                      Temperature = phys_oce.sub.2014$temp_deg_c,
                      Salinity = phys_oce.sub.2014$salinity,
                      Fluorometry = phys_oce.sub.2014$flurom_arbit,
                      "N Species" = phys_oce.sub.2014$nitrogen_species,
                      Silicate = phys_oce.sub.2014$SiOH4_mumol_l,
                      Phosphate = phys_oce.sub.2014$PO4_mumol_l,
                      Latitude = stat_names.2014$latitude,
                      Longitude = stat_names.2014$longitude))

corrdat.long <- cor(cbind(Depth = phys_oce.sub.long$depth,
                      Temperature = phys_oce.sub.long$temp_deg_c,
                      Salinity = phys_oce.sub.long$salinity,
                      Fluorometry = phys_oce.sub.long$flurom_arbit,
                      "N Species" = phys_oce.sub.long$nitrogen_species,
                      Silicate = phys_oce.sub.long$SiOH4_mumol_l,
                      Phosphate = phys_oce.sub.long$PO4_mumol_l,
                      Year = stat_names.long$year,
                      Latitude = stat_names.long$latitude,
                      Longitude = stat_names.long$longitude))

corrdat.2016 <- cor(cbind(Depth = phys_oce.sub.2016$depth,
                      Temperature = phys_oce.sub.2016$temp_deg_c,
                      Salinity = phys_oce.sub.2016$salinity,
                      Fluorometry = phys_oce.sub.2016$flurom_arbit,
                      "Ice Coverage" = phys_oce.sub.2016$icecover,
                      Latitude = stat_names.2016$latitude,
                      Longitude = stat_names.2016$longitude))


library(corrplot)

png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/2014_corr1.png")
corrplot(corrdat.2014, 
         method = "ellipse",
         type = "lower")
dev.off()
png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/2014_corr2.png")
corrplot(corrdat.2014, 
         method = "number",
         type = "lower")
dev.off()

png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/long_corr1.png")
corrplot(corrdat.long, 
         method = "ellipse",
         type = "lower")
dev.off()
png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/long_corr2.png")
corrplot(corrdat.long, 
         method = "number",
         type = "lower")
dev.off()

png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/2016_corr1.png")
corrplot(corrdat.2016, 
         method = "ellipse",
         type = "lower")
dev.off()
png("~/Studium/19SS/BA/ba_thesis_report/appendix_plots/2016_corr2.png")
corrplot(corrdat.2016, 
         method = "number",
         type = "lower")
dev.off()
