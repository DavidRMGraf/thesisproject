# script correlations physical data
rm(list=ls())
graphics.off()


phys_dat <- read.csv("physical_data_matrix.csv")[,c(1,12:22)]
str(phys_dat)
attach(phys_dat)
cor.test(depth_water_ctd, pressure_dbar)
# pressure kann raus
cor.test(depth_water_ctd, temp_deg_c)
# die Korrelation zw. Tiefe und Temperatur ist n.s.

cor.test(temp_deg_c, temp_pot_deg_c)
# temp pot kann raus

cor.test(conductivity_ms_cm, salinity)
# conductivity kann raus

cor.test(temp_deg_c, salinity)
# starke korrelation auch bei salinity und temperature

cor.test(flurom_arbit, temp_deg_c)
plot(temp_deg_c, flurom_arbit)
abline(lm(flurom_arbit~temp_deg_c))
hist(flurom_arbit)

range(flurom_arbit)
# fluorescence is supposed to be a proxy for ChlA -> depends on wavelength!

detach(phys_dat)

#' auf Basis hiervon würde ich Tiefe, salinity, temperature und flurometry behalten