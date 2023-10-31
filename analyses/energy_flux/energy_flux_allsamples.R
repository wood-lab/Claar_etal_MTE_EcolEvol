# Install necessary libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork)

# Load necessary data
# data <- read.csv("data/historical_data/stickle_para_data_2020-01-09.csv")
load("data/historical_data/historical_data.RData")

# Standardize weights from hundredths of a gram to grams
data$wt_schist <- data$wt_schist*0.01
data$weight <- data$weight*0.01
data$num_schist_1_2plus <- data$num_schist
data$num_schist_1_2plus[data$num_schist_1_2plus>1] <- "2+"

# Select only infected fish 
schist1 <- data %>% 
  filter(num_schist>0) %>%
  filter(weight != "NA") %>% 
  filter(num_schist != "NA")

# Calculate n by littoral/limnetic
schist1 %>% group_by(gear) %>% tally()

# Pull out weights of each 
schisto.weights <- schist1$wt_schist
stickle.weights <- (schist1$weight - schist1$wt_schist)
number.schisto <- schist1$num_schist
stickle.length <- schist1$length
sort(schisto.weights)

# Create function for flux equation; set default values
flux_schisto <- function(i=exp(17.17), # Brown et al 2001 inverts  
                         Mp, 
                         alphap=3/4, # theoretical value
                         E=0.63, # theoretical value
                         Temp=273.15+15,  # in Kelvin, add 273.15 to value in C
                         Np=1){ # number of parasites
  Fp = (i*Mp^(alphap)*exp(-E/(8.62*10^(-5)*Temp)))*Np
  Fp
}

flux_stickle <- function(i=exp(17.71), # lm3 on 7/11/21
                         Mp, 
                         alphap=1.039, # lm3 on 7/11/21
                         E=0.63, # theoretical value
                         Temp=273.15+15,  # in Kelvin, add 273.15 to value in C
                         Np=1){ # number of parasites
  Fp = (i*Mp^(alphap)*exp(-E/(8.62*10^(-5)*Temp)))*Np
  Fp
}

# Calculate instantaneous energy flux for each schisto
schisto.flux <- flux_schisto(Mp=schisto.weights) %>% as.data.frame()

# Calculate instantaneous energy flux for each stickle
stickle.flux <- flux_stickle(Mp=stickle.weights) %>% as.data.frame()

# Calculate the amount of the stickle's energy going to the schisto
energy.siphoned <- schisto.flux

# Calculate the percent of the stickle's energy going to the schisto
energy.siphoned.percent <- (schisto.flux/stickle.flux)*100

# Calulate proportion of worm weight to stickle weight
schisto.weight.percent <- schisto.weights/stickle.weights*100

# Calulate difference in proportions of worm weight & stickle weight
schisto.weight.diff <- schisto.weight.percent-energy.siphoned.percent

flux_df <- cbind(schisto.weights,stickle.weights,schisto.flux,stickle.flux,energy.siphoned,energy.siphoned.percent,schist1$num_schist_1_2plus,number.schisto,stickle.length,schisto.weight.percent,schisto.weight.diff)
colnames(flux_df) <- c("schisto.weights","stickle.weights","schisto.flux","stickle.flux","energy.siphoned","energy.siphoned.perc","num_schist_1_2plus","number.schisto","stickle.length","schisto.weight.perc","schisto.weight.diff")

flux_df_1only <- flux_df %>% 
  filter(number.schisto==1)

flux_df_2plus <- flux_df %>% 
  filter(number.schisto>1)

max_Esip <- max(flux_df$energy.siphoned)
min_Esip <- min(flux_df$energy.siphoned)
max_fishwt <- max(flux_df$stickle.weights)
min_fishwt <- min(flux_df$stickle.weights)
max_schistwt <- max(flux_df$schisto.weights)
min_schistwt <- 0.0003 # manually set so it matches iliamna minimum weight
min_Esipper <- min(flux_df$energy.siphoned.perc)
max_Esipper <- max(flux_df$energy.siphoned.perc)

# Use for ms
min_Esipper_1only <- min(flux_df_1only$energy.siphoned.perc)
max_Esipper_1only <- max(flux_df_1only$energy.siphoned.perc)
max_Esip_1only <- max(flux_df_1only$energy.siphoned)
min_Esip_1only <- min(flux_df_1only$energy.siphoned)
mean_Esip_1only <- mean(flux_df_1only$energy.siphoned)
sd_Esip_1only <- sd(flux_df_1only$energy.siphoned)
mean_Esip_2plus <- mean(flux_df_2plus$energy.siphoned)
sd_Esip_2plus <- sd(flux_df_2plus$energy.siphoned)
min_Esipper_2plus <- min(flux_df_2plus$energy.siphoned.perc)
max_Esipper_2plus <- max(flux_df_2plus$energy.siphoned.perc)
max_Esip_2plus <- max(flux_df_2plus$energy.siphoned)
min_Esip_2plus <- min(flux_df_2plus$energy.siphoned)



p1 <- ggplot(flux_df_1only)+ theme_classic()+
  theme(legend.position = c(0.85,0.65))+
  geom_point(mapping=aes(x=stickle.weights,y=energy.siphoned.perc,
                         color=schisto.weights),shape=16)+
  scale_x_continuous(name = "host weight (g)",
                     limits=c(0,max_fishwt)) +
  scale_y_continuous(name = "energy siphoned (%)",
                     limits=c(0,max_Esipper)) +
  scale_color_viridis(option="viridis",direction=-1, end=0.8,
                      limits=c(min_schistwt, max_schistwt))+
  labs(color="parasite
weight (g)")
p1


p2<- ggplot(flux_df_1only, aes(x = stickle.weights, y = schisto.weights, 
                         color = energy.siphoned)) + 
  theme_classic() + 
  theme(legend.position = c(0.9,0.65),
        legend.background = element_blank())+
  geom_point(shape=16) + 
  scale_x_continuous(name = "host weight (g)",
                     limits=c(0,max_fishwt)) +
  scale_y_continuous(name = "parasite weight (g)",
                     limits=c(0,max_schistwt)) +
  scale_color_viridis(option="magma",direction=-1, end=0.85,
                      limits=c(min_Esip, max_Esip))+  
  labs(color="energy
siphoned")
p2


p3 <- ggplot(flux_df_2plus)+ theme_classic()+
  theme(legend.position = "none")+
  geom_point(mapping=aes(x=stickle.weights,y=energy.siphoned.perc,
                         color=schisto.weights),shape=17)+
  scale_x_continuous(name = "host weight (g)",
                     limits=c(0,max_fishwt)) +
  scale_y_continuous(name = "energy siphoned (%)",
                     limits=c(0,max_Esipper)) +
  scale_color_viridis(option="viridis",direction=-1, end=0.8,
                      limits=c(min_schistwt, max_schistwt))+
  labs(color="parasite
weight (g)")
p3

p4<- ggplot(flux_df_2plus, aes(x = stickle.weights, 
                               y = schisto.weights, 
                               color = energy.siphoned)) + 
  theme_classic() + 
  theme(legend.position = "none")+
  geom_point(shape=17) + 
  scale_x_continuous(name = "host weight (g)",
                     limits=c(0,max_fishwt)) +
  scale_y_continuous(name = "parasite weight (g)",
                     limits=c(0,max_schistwt)) +
  scale_color_viridis(option="magma",direction=-1,end=0.85,
                      limits=c(min_Esip, max_Esip) )+  
  labs(color="energy
siphoned")
p4

# 
# 
# 
# jpeg(file="analyses/energy_flux/energy_flux_multipanel_allsamples.jpg",
#     width = 8, height= 8,units="in",res=300) # open jpg image
# (p1+p3)/(p2+p4)
# dev.off()

# jpeg(file="analyses/energy_flux/energy_flux_multipanel_allsamples_panel2.jpg", width = 4, height= 4,units="in",res=300) # open jpg image
# p3 + scale_x_continuous(name = "host weight (g)") +
#   scale_y_continuous(name = "energy siphoned (%)") +
#   labs(color="parasite
# weight (g)")
# dev.off()

p5 <- ggplot(flux_df_1only)+ theme_classic()+
  theme(legend.position = "none")+
  geom_point(mapping=aes(x=stickle.length,y=energy.siphoned.perc,
                         color=schisto.weights),shape=16)+  
  scale_y_continuous(name = "energy siphoned (%)",
                     limits=c(0,max_Esipper)) +
  scale_color_viridis(option="viridis",direction=-1, end=0.8,
                      limits=c(min_schistwt, max_schistwt))+
  labs(color="parasite
weight (g)")
p5

p6 <- ggplot(flux_df)+theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name = "energy siphoned (%)",
                     limits=c(0,max_Esipper)) +
  scale_x_continuous(name="parasite:host mass (%)")+
  geom_point(mapping=aes(x=schisto.weight.perc,
                         y=energy.siphoned.perc,
                         color=schisto.weights),shape=16)+ 
  geom_abline(slope=1, intercept=0,cex=2)+
  scale_color_viridis(option="viridis",direction=-1, end=0.8,
                      limits=c(min_schistwt, max_schistwt))+
  NULL
p6

p7 <-ggplot(flux_df)+theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="parasite:host mass (%) - energy siphoned %")+
  scale_x_continuous(name="parasite:host mass (%)")+
  geom_point(mapping=aes(x=schisto.weight.perc,
                         y=schisto.weight.diff,
                         color=energy.siphoned),shape=16)+ 
  geom_hline(yintercept=0,cex=2)+
  scale_color_viridis(option="magma",direction=-1,end=0.85,
                      limits=c(min_Esip, max_Esip) )+  
  NULL
p7

p6+p7

jpeg(file="analyses/energy_flux/energy_flux_multipanel_allsamples_pluswtratio.jpg",
     width = 8, height= 12,units="in",res=300) # open jpg image
(p1+p2)/(p3+p4)/(p6+p7)+plot_annotation(tag_levels = "A")
dev.off()

pdf(file="analyses/energy_flux/energy_flux_multipanel_allsamples_pluswtratio.pdf",useDingbats = FALSE,
     width = 8, height= 12) # open  image
(p1+p2)/(p3+p4)/(p6+p7)+plot_annotation(tag_levels = "A")
dev.off()

p8 <-ggplot(flux_df)+theme_classic()+
  theme(legend.position = "none")+
  geom_point(mapping=aes(x=stickle.flux,
                         y=schisto.flux,
                         color=energy.siphoned),shape=16)+ 
  geom_abline(slope=1, cex=2)+
  scale_color_viridis(option="magma",direction=-1,end=0.85)+  
  NULL
p8

max(flux_df_1only$energy.siphoned.perc)
min(flux_df_1only$energy.siphoned.perc)
mean(flux_df_1only$energy.siphoned.perc)
sd(flux_df_1only$energy.siphoned.perc)

max(flux_df_2plus$energy.siphoned.perc)
min(flux_df_2plus$energy.siphoned.perc)
mean(flux_df_2plus$energy.siphoned.perc)
sd(flux_df_2plus$energy.siphoned.perc)

max(flux_df_1only$energy.siphoned)
min(flux_df_1only$energy.siphoned)
mean(flux_df_1only$energy.siphoned)
sd(flux_df_1only$energy.siphoned)

max(flux_df_2plus$energy.siphoned)
min(flux_df_2plus$energy.siphoned)
mean(flux_df_2plus$energy.siphoned)
sd(flux_df_2plus$energy.siphoned)


max(flux_df_1only$schisto.weight.perc)
min(flux_df_1only$schisto.weight.perc)
mean(flux_df_1only$schisto.weight.perc)
sd(flux_df_1only$schisto.weight.perc)

max(flux_df_2plus$schisto.weight.perc)
min(flux_df_2plus$schisto.weight.perc)
mean(flux_df_2plus$schisto.weight.perc)
sd(flux_df_2plus$schisto.weight.perc)

#####
lm10 <- lm(schisto.flux~stickle.flux,data=flux_df)
summary(lm10)
lm10$coefficients["stickle.flux"]
mean(flux_df$energy.siphoned.perc)
sd(flux_df$energy.siphoned.perc)


flux_df2 <- flux_df %>% filter(schisto.flux<0.002)
lm12 <- lm(schisto.flux~stickle.flux,data=flux_df2)
summary(lm12)
