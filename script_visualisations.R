## vizualisation script for results
# source the analysis script to get the R environment:
rm(list=ls())
graphics.off()
source(file = "script_analysis_all_data.r")
library(gridExtra)

## 2016 -----------
dm.plotdata.2016 <- cbind(low.sequ.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)
pca.plotdata.2016 <- cbind(coords.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)

## minus_1
p1 <- ggplot(data = dm.plotdata.2016, aes(x = -depth, y = minus_1))+
  geom_point();p1
## minus_2
p2 <- ggplot(data = dm.plotdata.2016, aes(x = temp_deg_c-2, y = minus_2))+
  geom_point();p2

p3 <- ggplot(data = dm.plotdata.2016, aes(x = icecover, y = minus_2))+
  geom_point();p3
## minus_3 not conclusive
## minus_4 not conclusive either


fviz_eig(pca.sequ.2016, addlabels = TRUE, ylim = c(0, 50), main = "Scree plot - 2016 subset")

# Dim.1
p1 <- ggplot(pca.plotdata.2016, aes(x = depth, y = Dim.1))+
  geom_point();p1
## Dim.2/Dim.3 - temp vs. icecov
p2 <- ggplot(pca.plotdata.2016, aes(x = temp_deg_c-2, y = Dim.2))+
  geom_point();p2
p3 <- ggplot(pca.plotdata.2016, aes(x = temp_deg_c-2, y = Dim.3))+
  geom_point();p3

p4 <- ggplot(pca.plotdata.2016, aes(x = icecover, y = Dim.2))+
  geom_point();p4
p5 <- ggplot(pca.plotdata.2016, aes(x = icecover, y = Dim.3))+
  geom_point();p5


## 2014-----
graphics.off()
dm.plotdata.2014 <- cbind(low.sequ.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
pca.plotdata.2014 <- cbind(coords.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)

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

fviz_eig(pca.sequ.2014, addlabels = TRUE, ylim = c(0, 50), main = "Scree plot - 2014 subset")
## Dim.1
p1 <- ggplot(pca.plotdata.2014, aes(x = longitude, y = Dim.1))+
  geom_point();p1
p1 <- ggplot(pca.plotdata.2014, aes(x = temp_deg_c-2, y = Dim.1))+
  geom_point();p1

## long -------
dm.plotdata.long <- cbind(low.sequ.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
pca.plotdata.long <- cbind(coords.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)

## minus_1
p1 <- ggplot(data = dm.plotdata.long, aes(y = temp_deg_c-2, x = year))+
  geom_point();p1

p2 <- ggplot(data = dm.plotdata.long, aes(x = year, y = minus_1))+
  geom_point();p2

fviz_eig(pca.sequ.long, addlabels = TRUE, ylim = c(0, 50), main = "Scree plot - long subset")
## Dim.1
p1 <- ggplot(pca.plotdata.long, aes(x = temp_deg_c-2, y = Dim.1))+
  geom_point();p1
## Dim.2
p2 <- ggplot(pca.plotdata.long, aes(x = year, y = Dim.2))+
  geom_point();p2
## Dim.3
p3 <- ggplot(pca.plotdata.long, aes(x = SiOH4_mumol_l, y = Dim.3))+
  geom_point();p3


### clustering -------------
for(i in 1:ncol(phys_oce.sub.2014)){
  phys_oce.sub.2014[,i] <- (phys_oce.sub.2014[,i] - mean(phys_oce.sub.2014[,i]))/sd(phys_oce.sub.2014[,i])
}

dm.plotdata.2014$clust <- as.factor(kmeans(phys_oce.sub.2014, centers = 2, nstart = 25)$cluster)
pca.plotdata.2014$clust <- dm.plotdata.2014$clust

p1 <- ggplot(data = dm.plotdata.2014, aes(x = longitude, y = latitude, col = clust))+
  geom_point();p1

p2 <- ggplot(data = pca.plotdata.2014, aes(x = temp_deg_c-2, y = Dim.1, col = clust))+
  geom_point();p2

