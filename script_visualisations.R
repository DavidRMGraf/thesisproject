## vizualisation script for results
# source the analysis script to get the R environment:
rm(list=ls())
graphics.off()
source(file = "script_analysis_all_data.r")
library(gridExtra)
#library(Cairo)
library(stickylabeller)
## 2014-----

dm.plotdata.2014 <- cbind(low.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
pca.plotdata.2014 <- cbind(coords.2014, phys_oce.sub.2014, year = stat_names.2014$year, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## correl 2014
vars.2014 <- cbind(phys_oce.sub.2014, latitude = stat_names.2014$latitude, longitude = stat_names.2014$longitude)
dm.correl.2014 <- as.data.frame.matrix(matrix(nrow = ncol(vars.2014), ncol = 12))
colnames(dm.correl.2014) <- rep(c("cor. coeff", "R.squ.", "pval"), 4)
rownames(dm.correl.2014) <- colnames(vars.2014)
pca.correl.2014 <- dm.correl.2014

for(i in 1:ncol(vars.2014)){
  for(j in 1:4){
    dm.correl.2014[i, (3*j)-2] <- round(cor(vars.2014[,i], low.2014[,j]), 2)
    dm.correl.2014[i, (3*j)-1] <- round(summary(lm(low.2014[,j]~vars.2014[,i]))$adj.r.squared, 2)
    dm.correl.2014[i, 3*j] <- summary(lm(low.2014[,j]~vars.2014[,i]))$coefficients[2,4]
    pca.correl.2014[i, (3*j)-2] <- round(cor(vars.2014[,i], pca.plotdata.2014[,j]), 2)
    pca.correl.2014[i, (3*j)-1] <- round(summary(lm(pca.plotdata.2014[,j]~vars.2014[,i]))$adj.r.squared, 2)
    pca.correl.2014[i, 3*j] <- summary(lm(pca.plotdata.2014[,j]~vars.2014[,i]))$coefficients[2,4]
  }
}
pca.correl.2014 <- pca.correl.2014[,-(7:12)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DM results
## minus_1

minus_1.2014 <- data.frame(minus_1 = rep(dm.plotdata.2014$minus_1, 4), 
                           value = c(dm.plotdata.2014$salinity, dm.plotdata.2014$longitude, dm.plotdata.2014$temp_deg_c, dm.plotdata.2014$latitude),
                           var = factor(c(rep("Salinity, R² = 0.31", 42), rep("Longitude, R² = 0.21", 42),  rep("Temperature, R² = 0.20", 42), rep("Latitude, R² = 0.18", 42)), 
                                        levels = c("Salinity, R² = 0.31", "Longitude, R² = 0.21", "Temperature, R² = 0.20", "Latitude, R² = 0.18"))
)

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2014_dc1.png", width = 1200, height = 750)
ggplot(minus_1.2014, aes(value, minus_1))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Value",
       y = "Diffusion Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free",
             labeller = label_glue('{.l}) {var}'))
dev.off()

## minus_2
minus_2.2014 <- data.frame(minus_2 = rep(dm.plotdata.2014$minus_2, 2), 
                           value = c(dm.plotdata.2014$depth, dm.plotdata.2014$longitude),
                           var = factor(c(rep("Depth, R² = 0.21", 42),  rep("Longitude, R² = 0.20", 42)), 
                                        levels = c("Depth, R² = 0.21", "Longitude, R² = 0.20"))
)

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2014_dc2.png", width = 1140, height = 570)
ggplot(minus_2.2014, aes(value, minus_2))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Value",
       y = "Diffusion Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## minus_3
png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2014_dc3.png", width = 600, height = 450)
ggplot(dm.plotdata.2014, aes(x = nitrogen_species, minus_3))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Nitrogen Species, R² = 0.16",
       y = "Diffusion Component 3")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()

## minus_4
png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2014_dc4.png", width = 600, height = 450)
ggplot(dm.plotdata.2014, aes(x = salinity, minus_4))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Salinity, R² = 0.24",
       y = "Diffusion Component 3")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PCA results
## Dim.1
## only first two dims
comp_1.2014 <- data.frame(Dim.1 = rep(pca.plotdata.2014$Dim.1, 4), 
                           value = c(pca.plotdata.2014$salinity, pca.plotdata.2014$longitude, pca.plotdata.2014$temp_deg_c, pca.plotdata.2014$latitude),
                          var = factor(c(rep("Salinity, R² = 0.32", 42), rep("Longitude, R² = 0.32", 42),  rep("Temperature, R² = 0.29", 42), rep("Latitude, R² = 0.21", 42)), 
                                       levels = c("Salinity, R² = 0.32", "Longitude, R² = 0.32", "Temperature, R² = 0.29", "Latitude, R² = 0.21"))
)

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/2014_pc1.png", width = 1200, height = 750)
ggplot(comp_1.2014, aes(value, Dim.1))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Value",
       y = "Principal Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## Dim.2

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/2014_pc2.png", width = 600, height = 452)
p2 <- ggplot(pca.plotdata.2014, aes(x = depth, y = Dim.2))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Depth, R² = 0.28",
       y = "Principal Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8);p2
dev.off()

## long -------
dm.plotdata.long <- cbind(low.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)
pca.plotdata.long <- cbind(coords.long, phys_oce.sub.long, year = stat_names.long$year, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## correl long
vars.long <- cbind(phys_oce.sub.long, latitude = stat_names.long$latitude, longitude = stat_names.long$longitude, year = stat_names.long$year)
dm.correl.long <- as.data.frame.matrix(matrix(nrow = ncol(vars.long), ncol = 12))
colnames(dm.correl.long) <- rep(c("cor. coeff", "R.squ.", "pval"), 4)
rownames(dm.correl.long) <- colnames(vars.long)
pca.correl.long <- dm.correl.long

for(i in 1:ncol(vars.long)){
  for(j in 1:4){
    dm.correl.long[i, (3*j)-2] <- round(cor(vars.long[,i], low.long[,j]), 2)
    dm.correl.long[i, (3*j)-1] <- round(summary(lm(low.long[,j]~vars.long[,i]))$adj.r.squared, 2)
    dm.correl.long[i, 3*j] <- summary(lm(low.long[,j]~vars.long[,i]))$coefficients[2,4]
    pca.correl.long[i, (3*j)-2] <- round(cor(vars.long[,i], pca.plotdata.long[,j]), 2)
    pca.correl.long[i, (3*j)-1] <- round(summary(lm(pca.plotdata.long[,j]~vars.long[,i]))$adj.r.squared, 2)
    pca.correl.long[i, 3*j] <- summary(lm(pca.plotdata.long[,j]~vars.long[,i]))$coefficients[2,4]
  }
}
pca.correl.long <- pca.correl.long[,-(10:12)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PCA results
## first three dims
## Dim.1
# temperature (!) and salinity, longitude (+/-)
comp_1.long <- data.frame(Dim.1 = rep(pca.plotdata.long$Dim.1, 3), 
                          value = c(pca.plotdata.long$temp_deg_c, pca.plotdata.long$salinity, pca.plotdata.long$longitude),
                          var = factor(c(rep("Temperature, R² = 0.79", 46), rep("Salinity, R² = 0.47", 46), rep("Longitude, R² = 0.34", 46)), 
                                       levels = c("Temperature, R² = 0.79", "Salinity, R² = 0.47", "Longitude, R² = 0.34")))

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/long_pc1.png", width = 1200, height = 400)
ggplot(comp_1.long, aes(value, Dim.1))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Value",
       y = "Principal Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## Dim.2
## depth and year
comp_2.long <- data.frame(Dim.2 = rep(pca.plotdata.long$Dim.2, 2), 
                          value = c(pca.plotdata.long$depth, pca.plotdata.long$year),
                          var = factor(c(rep("Depth, R² = 0.30", 46), rep("Year, R² = 0.21", 46)), 
                                       levels = c("Depth, R² = 0.30", "Year, R² = 0.21")))

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/long_pc2.png", width = 900, height = 450)
ggplot(comp_2.long, aes(value, Dim.2))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Value",
       y = "Principal Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## Dim.3
comp_3.long <- data.frame(Dim.3 = rep(pca.plotdata.long$Dim.3, 4), 
                          value = c(pca.plotdata.long$salinity, pca.plotdata.long$SiOH4_mumol_l,  pca.plotdata.long$year,  pca.plotdata.long$longitude),
                          var = factor(c(rep("Salinity, R² = 0.14", 46), rep("Silicate, R² = 0.10", 46), rep("Year, R² = 0.09", 46), rep("Longitude, R² = 0.06", 46)), 
                                       levels = c("Salinity, R² = 0.14", "Silicate, R² = 0.10", "Year, R² = 0.09", "Longitude, R² = 0.06")))

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/long_pc3.png", width = 1200, height = 750)
ggplot(comp_3.long, aes(value, Dim.3))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Value",
       y = "Principal Component 3")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DM results
## minus_1
minus_1.long <- data.frame(minus_1 = rep(dm.plotdata.long$minus_1, 3), 
                           value = c(dm.plotdata.long$temp_deg_c, dm.plotdata.long$salinity, dm.plotdata.long$longitude),
                           var = factor(c(rep("Temperature, R² = 0.72", 46), rep("Salinity, R² = 0.37", 46),  rep("Longitude, R² = 0.26", 46)), 
                                        levels = c("Temperature, R² = 0.72", "Salinity, R² = 0.37", "Longitude, R² = 0.26"))
)

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/long_dc1.png", width = 1200, height = 450)
ggplot(minus_1.long, aes(value, minus_1))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Value",
       y = "Diffusion Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## minus_2

minus_2.long <- data.frame(minus_2 = rep(dm.plotdata.long$minus_2, 3), 
                           value = c(dm.plotdata.long$depth, dm.plotdata.long$longitude, dm.plotdata.long$SiOH4_mumol_l),
                           var = factor(c(rep("Depth, R² = 0.14", 46), rep("Silicate, R² = 0.07", 46), rep("Longitude, R² = 0.07", 46)), 
                                        levels = c("Depth, R² = 0.14", "Silicate, R² = 0.07", "Longitude, R² = 0.07"))
)

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/long_dc2.png", width = 1200, height = 450)
ggplot(minus_2.long, aes(value, minus_2))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Value",
       y = "Diffusion Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## minus_3

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/long_dc3.png", width = 600, height = 450)
ggplot(dm.plotdata.long, aes(year, minus_2))+
  geom_point(size = 3)+
  labs(x = "Year, R² = 0.31",
       y = "Diffusion Component 3")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()

## minus_4

## 2016 -----------
dm.plotdata.2016 <- cbind(low.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)
pca.plotdata.2016 <- cbind(coords.2016, phys_oce.sub.2016, year = stat_names.2016$year, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## correl
vars.2016 <- cbind(phys_oce.sub.2016, latitude = stat_names.2016$latitude, longitude = stat_names.2016$longitude)
dm.correl.2016 <- as.data.frame.matrix(matrix(nrow = ncol(vars.2016), ncol = 12))
colnames(dm.correl.2016) <- rep(c("cor. coeff", "R.squ.", "pval"), 4)
rownames(dm.correl.2016) <- colnames(vars.2016)
pca.correl.2016 <- dm.correl.2016

for(i in 1:ncol(vars.2016)){
  for(j in 1:4){
    dm.correl.2016[i, (3*j)-2] <- round(cor(vars.2016[,i], low.2016[,j]), 2)
    dm.correl.2016[i, (3*j)-1] <- round(summary(lm(low.2016[,j]~vars.2016[,i]))$adj.r.squared, 2)
    dm.correl.2016[i, 3*j] <- summary(lm(low.2016[,j]~vars.2016[,i]))$coefficients[2,4]
    pca.correl.2016[i, (3*j)-2] <- round(cor(vars.2016[,i], pca.plotdata.2016[,j]), 2)
    pca.correl.2016[i, (3*j)-1] <- round(summary(lm(pca.plotdata.2016[,j]~vars.2016[,i]))$adj.r.squared, 2)
    pca.correl.2016[i, 3*j] <- summary(lm(pca.plotdata.2016[,j]~vars.2016[,i]))$coefficients[2,4]
  }
}
pca.correl.2016 <- pca.correl.2016[,-(7:12)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DM results
## minus_1

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2016_dc1.png", width = 600, height = 450)
ggplot(dm.plotdata.2016, aes(depth, minus_1))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Depth, R² = 0.89",
       y = "Diffusion Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()

## minus_2

minus_2.2016 <- data.frame(minus_2 = rep(dm.plotdata.2016$minus_2, 4), 
                           value = c(dm.plotdata.2016$temp_deg_c, dm.plotdata.2016$longitude, dm.plotdata.2016$icecover, dm.plotdata.2016$salinity),
                           var = factor(c(rep("Temperature, R² = 0.82", 64), rep("Longitude, R² = 0.72", 64), rep("Ice Coverage, R² = 0.68", 64), rep("Salinity, R² = 0.31", 64)), 
                                        levels = c("Temperature, R² = 0.82", "Longitude, R² = 0.72", "Ice Coverage, R² = 0.68", "Salinity, R² = 0.31"))
)

png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2016_dc2.png", width = 1200, height = 750)
ggplot(minus_2.2016, aes(value, minus_2))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Value",
       y = "Diffusion Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()

## minus_3 not conclusive
## minus_4 
png("~/Studium/19SS/BA/ba_thesis_report/dm_plots/2016_dc4.png", width = 600, height = 450)
ggplot(dm.plotdata.2016, aes(flurom_arbit, minus_4))+
  geom_point(size = 3, pch = 16)+
  labs(x = "Fluorometry, R² = 0.23",
       y = "Diffusion Component 4")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PCA results
## only first two dims!
# Dim.1

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/2016_pc1.png", width = 600, height = 450)
ggplot(pca.plotdata.2016, aes(depth, Dim.1))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Depth, R² = 0.78",
       y = "Principal Component 1")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)
dev.off()

# p1 <- ggplot(pca.plotdata.2016, aes(x = depth, y = Dim.1))+
#   geom_point();p1
# summary(lm(pca.plotdata.2016$Dim.1~pca.plotdata.2016$depth))

## Dim.2/Dim.3 - temp vs. icecov

comp_2.2016 <- data.frame(Dim.2 = rep(pca.plotdata.2016$Dim.2, 4), 
                          value = c(pca.plotdata.2016$temp_deg_c, pca.plotdata.2016$longitude,  pca.plotdata.2016$icecover,  pca.plotdata.2016$salinity),
                          var = factor(c(rep("Temperature, R² = 0.89", 64), rep("Longitude, R² = 0.60", 64), rep("Ice Coverage, R² = 0.58", 64), rep("Salinity, R² = 0.53", 64)), 
                                       levels = c("Temperature, R² = 0.89", "Longitude, R² = 0.60", "Ice Coverage, R² = 0.58", "Salinity, R² = 0.53")))

png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/2016_pc2.png", width = 1200, height = 750)
ggplot(comp_2.2016, aes(value, Dim.2))+
  geom_point(size = 3, pch = 17)+
  labs(x = "Value",
       y = "Principal Component 2")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text =  element_text(size = 20))+
  geom_smooth(method = "lm", se = F, col = "black", lwd = .8)+
  facet_wrap(~var, scales = "free", labeller = label_glue('{.l}) {var}'))
dev.off()


# screeplots -----------

scree_2014 <- fviz_eig(pca.2014, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "2014")+
  geom_vline(xintercept = 2.5, col = "black", lty = "dashed", lwd = .7)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        title = element_text(size = 20))

scree_2016 <-   fviz_eig(pca.2016, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "2016")+
  geom_vline(xintercept = 2.5, col = "black", lty = "dashed", lwd = .7)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        title = element_text(size = 20))

scree_long <- fviz_eig(pca.long, geom = "bar", addlabels = TRUE, ggtheme = theme_bw(), ylim = c(0, 40), main = "long")+
  geom_vline(xintercept = 3.5, col = "black", lty = "dashed", lwd = .7)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        title = element_text(size = 20))


png("~/Studium/19SS/BA/ba_thesis_report/pca_plots/screeplots.png", width = 889, height = 679)
grid.arrange(scree_2014, scree_long, scree_2016, nrow = 1,
             left = textGrob("Percentage of explained Variance", gp=gpar(fontsize=18), rot = 90),
             bottom = textGrob("Dimensions", gp=gpar(fontsize=18)))
dev.off()

### clustering -------------
# for(i in 1:ncol(phys_oce.sub.2014)){
#   phys_oce.sub.2014[,i] <- (phys_oce.sub.2014[,i] - mean(phys_oce.sub.2014[,i]))/sd(phys_oce.sub.2014[,i])
# }
# 
# dm.plotdata.2014$clust <- as.factor(kmeans(phys_oce.sub.2014, centers = 3, nstart = 25)$cluster)
# pca.plotdata.2014$clust <- dm.plotdata.2014$clust
# 
# p1 <- ggplot(data = dm.plotdata.2014, aes(x = longitude, y = -depth, col = clust))+
#   geom_point();p1
# 
# p2 <- ggplot(data = pca.plotdata.2014, aes(x = temp_deg_c, y = Dim.1, col = clust))+
#   geom_point();p2

