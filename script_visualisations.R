## vizualisation script for results
# source the analysis script to get the R environment:
rm(list=ls())
graphics.off()
source(file = "script_analysis_all_data.r")
library(gridExtra)
## graphs -----------------
dm.plotdata.2014 <- cbind(low.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
dm.plotdata.long <- cbind(low.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
dm.plotdata.2016 <- cbind(low.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)

dm.plotdata.sequ.2014 <- cbind(low.sequ.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
dm.plotdata.sequ.long <- cbind(low.sequ.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
dm.plotdata.sequ.2016 <- cbind(low.sequ.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)


## 2016 -----------
## minus_1
p1 <- ggplot(data = dm.plotdata.2016, aes(x = -depth, y = minus_1))+
  geom_point();p1
p1 <- ggplot(data = dm.plotdata.sequ.2016, aes(x = -depth, y = minus_1))+
  geom_point();p1
## minus_2
p2 <- ggplot(data = dm.plotdata.sequ.2016, aes(x = temp_deg_c-2, y = minus_2))+
  geom_point();p2
p3 <- ggplot(data = dm.plotdata.sequ.2016, aes(x = icecover, y = minus_2))+
  geom_point();p3

## 2014-----
graphics.off()
## minus_1
p1 <- ggplot(data = dm.plotdata.2014, aes(x = longitude, y = minus_1))+
  geom_point();p1

p2 <- ggplot(data = dm.plotdata.2014, aes(x = latitude, y = minus_1))+
  geom_point();p2

p3 <- ggplot(data = dm.plotdata.2014, aes(x = temp_deg_c-2, y = minus_1))+
  geom_point();p3

p4 <- ggplot(data = dm.plotdata.2014, aes(x = salinity, y = minus_1))+
  geom_point();p4

p5 <- ggplot(data = dm.plotdata.2014, aes(x = flurom_arbit, y = minus_1))+
  geom_point();p5

p6 <- ggplot(data = dm.plotdata.2014, aes(x = nitrogen_species, y = minus_1))+
  geom_point();p6

## minus_2
p7 <- ggplot(data = dm.plotdata.2014, aes(x = flurom_arbit, y = minus_2))+
  geom_point();p7

## minus_3
# p8 <- ggplot(data = dm.plotdata.2014, aes(x = salinity, y = minus_3))+
#   geom_point();p8
# p9 <- ggplot(data = dm.plotdata.2014, aes(x = nitrogen_species, y = minus_3))+
#   geom_point();p9

## long -------
## minus_1
p1 <- ggplot(data = dm.plotdata.long, aes(y = temp_deg_c-2, x = year, color = -depth))+
  geom_point();p1

p2 <- ggplot(data = dm.plotdata.long, aes(x = temp_deg_c, y = minus_1))+
  geom_point();p2

p3 <- ggplot(data = dm.plotdata.long, aes(x = year, y = minus_1))+
  geom_point();p3

## PCA ---------------------
pca.plotdata.2014 <- cbind(coords.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
pca.plotdata.long <- cbind(coords.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
pca.plotdata.2016 <- cbind(coords.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)

## 2016
p1 <- ggplot(pca.plotdata.2016, aes(x = depth, y = Dim.1))+
  geom_point();p1

p2 <- ggplot(pca.plotdata.2016, aes(x = temp_deg_c-2, y = Dim.2))+
  geom_point();p2
p3 <- ggplot(pca.plotdata.2016, aes(x = temp_deg_c-2, y = Dim.3))+
  geom_point();p3

p4 <- ggplot(pca.plotdata.2016, aes(x = icecover, y = Dim.2))+
  geom_point();p4
p5 <- ggplot(pca.plotdata.2016, aes(x = icecover, y = Dim.3))+
  geom_point();p5

##2014
p1 <- ggplot(pca.plotdata.2014, aes(x = longitude, y = Dim.1))+
  geom_point();p1
p1 <- ggplot(pca.plotdata.2014, aes(x = temp_deg_c-2, y = Dim.1))+
  geom_point();p1

### clustering -------------
dm.plotdata.2014$clust <- as.factor(kmeans(phys_oce.sub.2014, centers = 3, nstart = 25)$cluster)

p1 <- ggplot(data = dm.plotdata.2014, aes(x = temp_deg_c-2, y = minus_1, col = clust))+
  geom_point();p1

