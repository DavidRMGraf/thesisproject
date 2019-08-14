## vizualisation script for results
# source the analysis script to get the R environment:
rm(list=ls())
graphics.off()
source(file = "script_analysis_all_data.r")
library(gridExtra)

## 2016 -----------
dm.plotdata.2016 <- cbind(low.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)
pca.plotdata.2016 <- cbind(coords.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)

## minus_1
p1 <- ggplot(data = dm.plotdata.2016, aes(x = depth, y = minus_1))+
  geom_point();p1
summary(lm(dm.plotdata.2016$minus_1~dm.plotdata.2016$depth))

## minus_2
p2 <- ggplot(data = dm.plotdata.2016, aes(x = temp_deg_c, y = minus_2))+
  geom_point();p2
summary(lm(dm.plotdata.2016$minus_2~dm.plotdata.2016$temp_deg_c))

p3 <- ggplot(data = dm.plotdata.2016, aes(x = icecover, y = minus_2))+
  geom_point();p3
summary(lm(dm.plotdata.2016$minus_2~dm.plotdata.2016$icecover))

#' as p2 has a higher adj.R^2, it is concluded that temperature is likely the right
#' controlling variable, rather than icecover -> temp!

## minus_3 not conclusive
## minus_4 not conclusive either

## only first two dims!
# Dim.1
p1 <- ggplot(pca.plotdata.2016, aes(x = depth, y = Dim.1))+
  geom_point();p1
summary(lm(pca.plotdata.2016$Dim.1~pca.plotdata.2016$depth))

## Dim.2/Dim.3 - temp vs. icecov
p2 <- ggplot(pca.plotdata.2016, aes(x = temp_deg_c, y = Dim.2))+
  geom_point();p2
summary(lm(pca.plotdata.2016$Dim.2~pca.plotdata.2016$temp_deg_c))
p4 <- ggplot(pca.plotdata.2016, aes(x = icecover, y = Dim.2))+
  geom_point();p4
summary(lm(pca.plotdata.2016$Dim.2~pca.plotdata.2016$icecover))
#' higher adj. R^2 of the lm for temp. here as well -> temp!

## 2014-----

dm.plotdata.2014 <- cbind(low.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
pca.plotdata.2014 <- cbind(coords.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)

## minus_1
p1 <- ggplot(data = dm.plotdata.2014, aes(x = longitude, y = minus_1))+
  geom_point();p1
summary(lm(dm.plotdata.2014$minus_1~dm.plotdata.2014$longitude))

p2 <- ggplot(data = dm.plotdata.2014, aes(x = latitude, y = minus_1))+
  geom_point();p2
summary(lm(dm.plotdata.2014$minus_1~dm.plotdata.2014$latitude))

p3 <- ggplot(data = dm.plotdata.2014, aes(x = temp_deg_c, y = minus_1))+
  geom_point();p3
summary(lm(dm.plotdata.2014$minus_1~dm.plotdata.2014$temp_deg_c))

p4 <- ggplot(data = dm.plotdata.2014, aes(x = salinity, y = minus_1))+
  geom_point();p4
summary(lm(dm.plotdata.2014$minus_1~dm.plotdata.2014$salinity))

## minus_2
p5 <- ggplot(data = dm.plotdata.2014, aes(x = depth, y = minus_2))+
  geom_point();p5
summary(lm(dm.plotdata.2014$minus_2~dm.plotdata.2014$depth))

## only first two dims
## Dim.1 -> longitude!
p1 <- ggplot(pca.plotdata.2014, aes(x = longitude, y = Dim.1))+
  geom_point();p1
summary(lm(pca.plotdata.2014$Dim.1~pca.plotdata.2014$longitude))

p1 <- ggplot(pca.plotdata.2014, aes(x = temp_deg_c, y = Dim.1))+
  geom_point();p1
summary(lm(pca.plotdata.2014$Dim.1~pca.plotdata.2014$temp_deg_c))
#' -> longitude has higher mult. R^2 than temp! -> longitude

## Dim.2 -> depth!
p2 <- ggplot(pca.plotdata.2014, aes(x = depth, y = Dim.2))+
  geom_point();p2
summary(lm(pca.plotdata.2014$Dim.2~pca.plotdata.2014$depth))
#' depth!

## long -------
dm.plotdata.long <- cbind(low.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
pca.plotdata.long <- cbind(coords.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)

## minus_1
p1 <- ggplot(data = dm.plotdata.long, aes(y = temp_deg_c-2, x = year))+
  geom_point();p1

p2 <- ggplot(data = dm.plotdata.long, aes(x = year, y = minus_1))+
  geom_point();p2

scree_long <- fviz_eig(pca.long, addlabels = TRUE, ylim = c(0, 40), main = "long")

## first three dims
## Dim.1
p1 <- ggplot(pca.plotdata.long, aes(x = temp_deg_c-2, y = Dim.1))+
  geom_point();p1
## Dim.2
p2 <- ggplot(pca.plotdata.long, aes(x = year, y = Dim.2))+
  geom_point();p2
## Dim.3
p3 <- ggplot(pca.plotdata.long, aes(x = SiOH4_mumol_l, y = Dim.3))+
  geom_point();p3


scree_2014 <- fviz_eig(pca.2014, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "2014")+
  geom_vline(xintercept = 2.5, col = "black", lty = "dotted")+
  theme(axis.title = element_blank())

scree_2016 <- fviz_eig(pca.2016, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "2016")+
  geom_vline(xintercept = 2.5, col = "black", lty = "dotted")+
  theme(axis.title = element_blank())

scree_long <- fviz_eig(pca.long, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "long")+
  geom_vline(xintercept = 3.5, col = "black", lty = "dotted")+
  theme(axis.title = element_blank())

grid.arrange(scree_2016, scree_2014, scree_long, nrow = 1,
             bottom = "Dimensions",
             left = "Percentage of explained variance")

### clustering -------------
for(i in 1:ncol(phys_oce.sub.2014)){
  phys_oce.sub.2014[,i] <- (phys_oce.sub.2014[,i] - mean(phys_oce.sub.2014[,i]))/sd(phys_oce.sub.2014[,i])
}

dm.plotdata.2014$clust <- as.factor(kmeans(phys_oce.sub.2014, centers = 3, nstart = 25)$cluster)
pca.plotdata.2014$clust <- dm.plotdata.2014$clust

p1 <- ggplot(data = dm.plotdata.2014, aes(x = longitude, y = -depth, col = clust))+
  geom_point();p1

p2 <- ggplot(data = pca.plotdata.2014, aes(x = temp_deg_c, y = Dim.1, col = clust))+
  geom_point();p2

