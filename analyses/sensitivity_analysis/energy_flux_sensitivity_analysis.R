# Load necessary packages
library(ggplot2)
library(gridExtra)
library(tidyverse)
#library(gapminder)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(emdbook)


#####Description of parameters 
# Fp = Ip*Np
# Fp is the energy flux
# Ip is the mean individual metabolic rate
# Np is the infrapopulation size

# Fp is the energy flux
## what we're calculating for
# i is the normalization constant 
## varies for organisms of different physiological types
# Mp is the body mass of the parasite
# alphap is the scaling exponent
## often 3/4
# E is the activation energy for enzymatic reactions 
## average 0.63 eV for oxidative respiration
# k is Boltzmann’s constant
## k = (8.62*10^(-5) eV*K^(-1))
# T is temperature in Kelvin
# Np is infrapopulation size (i.e., # of parasites)


#####

##### Set values for each parameter 
# Set values for each parameter 
# Set values for each flux parameter
C_best_schisto <- 17.17 # Brown 2001 invertebrates
C_best_stickle <- 17.61989 # lm4 on 12/15
C_min_schisto <- 15
C_max_schisto <- 19
i_min_schisto <- exp(C_min_schisto)
i_max_schisto <- exp(C_max_schisto)
i_best_schisto <- exp(C_best_schisto)
alphap_min <- 0.3 # minimum alpha value updated 11/11/20
alphap_max <- 1.31 # maximum i value found in literature as of 2/8
alphap_std <- 0.75 # standard value used 
E_min <- 0.2 # minimum E value found in literature as of 2/8
E_max <- 1.2 # maximum E value found in literature as of 2/8
E_std <- 0.65
Temp_min <- 273.15 + 0 # minimum lake temperature
Temp_max <- 273.15 + 25 # maximum lake temperature - rounded up 
Temp_std <- 273.15 + 15 # temperature used in respirometry
Mp_min <- 0.01 # mimimum parasite weight (single infection)
Mp_max <- 0.57 # maximum parasite weight (single infection)
Mp_mean <- 0.22 # mean parasite weight (single infection)
Temp_r <- seq(Temp_min,Temp_max,0.1) # Create range of temps from min to max
i_r_schisto <- lseq(i_min_schisto,i_max_schisto,100)
Mp_r <- seq(Mp_min,Mp_max,0.01)
alphap_r <- seq(alphap_min,alphap_max,0.01)
E_r <- seq(E_min,E_max,0.01)


# Create function for flux equation; set default values
flux <- function(i=i_best_schisto, # 
                 Mp=0.22, # mean schisto weight, in grams
                 alphap=3/4, # theoretical value
                 E=0.65, # theoretical value
                 # k=8.62*10^(-5), # Boltzmann's constant (doesn't change)
                 Temp=273.15,  # in Kelvin, add 273.15 to value in C
                 Np=1){ # number of parasites
  Fp = (i*Mp^(alphap)*exp(-E/(8.62*10^(-5)*Temp)))*Np
  Fp
}

# Set plotting parameters for all plots
myPalette <- viridis_pal(alpha = 1, begin = 0.25, end = 0.9,
                         direction = -1, option = "inferno")
sc <- scale_colour_gradientn(colours = myPalette(100), 
                             oob=scales::squish)


##### Temperature #####
## Calculate flux sensitivity
# Create empty data frame
plot_flux_temp <- function(temps=Temp_r,
                           base_temp=Temp_std,
                           alphap=alphap_std,
                           E=E_std,i=i_best_schisto,
                           Mp=Mp_mean,
                           xlab=c("temperature (C)"),
                           plot_or_slope="plot",
                           color_or_black="color"){
  flux_temp <- data.frame(Temp=Temp_r,flux=NA,slope=NA)
  j=2 # Initialize counter, j

flux_temp$flux[1] <- flux(Temp = temps[1],E=E,alphap=alphap,i=i,Mp=Mp)
  # Make for loop to compute flux equation for each temperature value
for(temp in temps[2:length(temps)]){ # for every temperature in the sequence
    flux_temp$flux[j] <- flux(Temp = temp,
                              E=E,
                              alphap=alphap,
                              i=i,
                              Mp=Mp) # run the flux function
    flux_temp$slope[j] <- (flux_temp$flux[j-1] - flux_temp$flux[j])/(flux_temp$Temp[j-1]-flux_temp$Temp[j])
    j <- j+1 # move to the next temperature
    }
  ## Plot temperature vs. flux difference from mean
  # Calculate flux at mean temperature
  flux_temp_comp <- flux(Temp=base_temp,
                         E=E,
                         alphap=alphap,
                         i=i,
                         Mp=Mp)
  # Calculate the differrence between mean flux and all fluxes within range
  flux_temp$flux_diff <- (flux_temp$flux-flux_temp_comp)/(flux_temp_comp)*100
  # Create new column with temperature in C for readability
  flux_temp$Temp_C <- flux_temp$Temp-273.15
  # Make plot
  if(color_or_black=="color"){
      plot_flux_temp_diff <- ggplot(flux_temp,
                                aes(x=Temp_C,
                                    y=flux_diff,
                                    color=abs(slope)))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      geom_hline(yintercept = 0,color='gray')+ 
      ylim(c(-150,200))+ 
      theme(axis.title.y = element_blank(),
          legend.position = c(0.8,0.8)) + 
      xlab(xlab)+
      sc
  }
  else if(color_or_black=="black") {
    plot_flux_temp_diff <- ggplot(flux_temp,
                                  aes(x=Temp_C,
                                      y=flux_diff))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      geom_hline(yintercept = 0,color='gray')+ 
      ylim(c(-150,200))+ 
      theme(axis.title.y = element_blank(),
            legend.position = c(0.8,0.8)) + 
      xlab(xlab)
  } else {
    print("Invalid entry - only accepts 'color' or 'black'")
  }
  
  df <- data.frame(min=min(flux_temp$slope,na.rm=TRUE),
                   max=max(flux_temp$slope,na.rm=TRUE))
  rownames(df) <- "temp"
  
  # Output
  if (plot_or_slope == "plot"){
    plot_flux_temp_diff
  } else if (plot_or_slope == "slope") {
    print(df)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
  
}

plot_flux_temp(temps=Temp_r,base_temp=Temp_std,
               plot_or_slope = "plot",
               color_or_black = "black")

##### i #####
plot_flux_i <- function(base_i=i_best_schisto,E=E_std,
                        alphap=alphap_std,Mp=Mp_mean,
                        Temp=Temp_std,
                        xlab=c("i (constant)"),
                        color_or_black="color",
                        plot_or_slope="plot"){
  flux_i <- data.frame(i=i_r_schisto,flux=NA,slope=NA)
  j=2
  
  flux_i$flux[1] <- flux(i=i_r_schisto[1],E=E,
                         alphap=alphap,Temp=Temp_std,
                         Mp=Mp)
  
  for(each_i in i_r_schisto[2:length(i_r_schisto)]){
    flux_i$flux[j] <- flux(i=i_r_schisto[j],
                           E=E,
                           alphap=alphap,
                           Temp = Temp_std,
                           Mp=Mp)
    flux_i$slope[j] <- (flux_i$flux[j-1] - flux_i$flux[j])/(flux_i$i[j-1]-flux_i$i[j])
    j <- j+1
  }

  flux_i_comp <- flux(i=base_i,E=E,alphap=alphap,Temp = Temp,Mp=Mp)
  flux_i$flux_diff <- (flux_i$flux-flux_i_comp)/(flux_i_comp)*100
  
  if(color_or_black=="color"){
    plot_flux_i_diff <- ggplot(flux_i,aes(x=i_r_schisto,
                                        y=flux_diff,
                                        color=abs(slope)))+
    theme_classic()+
    coord_cartesian(expand =FALSE)+
    geom_line()+
    geom_hline(yintercept = 0,color='gray')+ 
    ylim(c(-150,200))+
    theme(axis.title.y = element_blank(),
          legend.position = c(0.8,0.8)) + 
    xlab(xlab)+
    sc
  } else if(color_or_black=="black") {
    plot_flux_i_diff <- ggplot(flux_i,aes(x=i_r_schisto,
                                          y=flux_diff))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      geom_hline(yintercept = 0,color='gray')+ 
      ylim(c(-150,200))+
      theme(axis.title.y = element_blank(),
            legend.position = c(0.8,0.8)) + 
      xlab(xlab)
  } else {
  print("Invalid entry - only accepts 'color' or 'black'")
  }
  
  df <- data.frame(min=min(flux_i$slope,na.rm=TRUE),
                   max=max(flux_i$slope,na.rm=TRUE))
  rownames(df) <- "i"
  
  # Output
  if (plot_or_slope == "plot"){
    plot_flux_i_diff
  } else if (plot_or_slope == "slope") {
    print(df)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
}
  
plot_flux_i(plot_or_slope = "plot")
plot_flux_i(plot_or_slope = "plot",color_or_black = "black")


##### Mp #####
plot_flux_Mp <- function(base_Mp=Mp_mean,
                         alphap=alphap_std,E=E_std,
                         i=i_best_schisto,
                         Temp=Temp_std,
                         xlab=c("parasite mass (g)"),
                         plot_or_slope="plot",
                         color_or_black="color"){
  flux_Mp <- data.frame(Mp=Mp_r,flux=NA,slope=NA)
  j=2
  flux_Mp$flux[1] <- flux(Mp = Mp_r[1],E=E,alphap=alphap,i=i,Temp=Temp_std)
  for(mass in Mp_r[2:length(Mp_r)]){
    flux_Mp$flux[j] <- flux(Mp=mass,alphap=alphap,
                            E=E,i=i,Temp = Temp)
    
    flux_Mp$slope[j] <- (flux_Mp$flux[j-1] - flux_Mp$flux[j])/(flux_Mp$Mp[j-1]-flux_Mp$Mp[j])
    
    j <- j+1
  }

  flux_Mp_comp <- flux(Mp=base_Mp,alphap=alphap,
                       E=E,i=i,Temp = Temp)
  flux_Mp$flux_diff <- (flux_Mp$flux-flux_Mp_comp)/(flux_Mp_comp)*100
  if(color_or_black=="color"){
    plot_flux_Mp_diff <- ggplot(flux_Mp,
                              aes(x=Mp_r,y=flux_diff,
                                  color=abs(slope)))+
    theme_classic()+
    coord_cartesian(expand =FALSE)+
    geom_line()+
    geom_hline(yintercept = 0,color='gray')+
    ylim(c(-150,200))+
    theme(axis.title.y = element_blank(),
          legend.position = c(0.8,0.8)) + 
    xlab(xlab)+
    sc
  } else if(color_or_black=="black") {
    plot_flux_Mp_diff <- ggplot(flux_Mp,
                                aes(x=Mp_r,y=flux_diff))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      geom_hline(yintercept = 0,color='gray')+
      ylim(c(-150,200))+
      theme(axis.title.y = element_blank(),
            legend.position = c(0.8,0.8)) + 
      xlab(xlab)
  } else {
    print("Invalid entry - only accepts 'color' or 'black'")
  }
  df <- data.frame(min=min(flux_Mp$slope,na.rm=TRUE),
                   max=max(flux_Mp$slope,na.rm=TRUE))
  rownames(df) <- "Mp"
  
  # Output
  if (plot_or_slope == "plot"){
    plot_flux_Mp_diff
  } else if (plot_or_slope == "slope") {
    print(df)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
  
  
}
plot_flux_Mp(plot_or_slope = "plot")
plot_flux_Mp(plot_or_slope = "plot",color_or_black = "black")

##### alpha #####
plot_flux_alphap <- function(base_alphap=alphap_std,
                             E=E_std,Mp=Mp_mean,
                             i=i_best_schisto,
                             Temp=Temp_std,
                             xlab=c("alpha (scaling exp)"),
                             plot_or_slope="plot",
                             color_or_black="color"){
  flux_alphap <- data.frame(alphap=alphap_r,
                            flux=NA,slope=NA)
  j=2
  
  flux_alphap$flux[1] <- flux(Temp = Temp_std,E=E,
                            alphap=alphap_r[1],i=i,Mp=Mp)
  
  for(alpha in alphap_r[2:length(alphap_r)]){
    flux_alphap$flux[j] <- flux(alphap=alpha,Mp=Mp,
                                i=i,Temp = Temp)
    flux_alphap$slope[j] <- (flux_alphap$flux[j-1] - flux_alphap$flux[j])/(flux_alphap$alphap[j-1]-flux_alphap$alphap[j])
    j <- j+1
  }

  flux_alphap_comp <- flux(alphap=base_alphap,Mp=Mp,i=i,Temp = Temp)
  flux_alphap$flux_diff <- (flux_alphap$flux-flux_alphap_comp)/(flux_alphap_comp)*100
  
  if(color_or_black=="color"){
  plot_flux_alphap_diff <- ggplot(flux_alphap,
                                  aes(x=alphap_r,
                                      y=flux_diff,
                                      color=abs(slope)))+
    theme_classic()+
    coord_cartesian(expand =FALSE)+
    geom_line()+
    theme(legend.position = c(0.8,0.8))+
    geom_hline(yintercept = 0,color='gray')+
    ylim(c(-150,200)) + 
    ylab("Flux (% diff)") + 
    xlab(xlab)+
    sc
  } else if(color_or_black=="black") {
    plot_flux_alphap_diff <- ggplot(flux_alphap,
                                    aes(x=alphap_r,
                                        y=flux_diff))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      theme(legend.position = c(0.8,0.8))+
      geom_hline(yintercept = 0,color='gray')+
      ylim(c(-150,200)) + 
      ylab("Flux (% diff)") + 
      xlab(xlab)    
  } else {
    print("Invalid entry - only accepts 'color' or 'black'")
  }

  df <- data.frame(min=min(flux_alphap$slope,na.rm=TRUE),
                   max=max(flux_alphap$slope,na.rm=TRUE))
  rownames(df) <- "alphap"
  
  # Output
  if (plot_or_slope == "plot"){
    plot_flux_alphap_diff
  } else if (plot_or_slope == "slope") {
    print(df)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
  
  }
plot_flux_alphap(plot_or_slope = "plot")
plot_flux_alphap(Mp=Mp_min,xlab="")
plot_flux_alphap(plot_or_slope = "plot",color_or_black = "black")


##### E #####
plot_flux_E <- function(base_E=E_std,Temp=Temp_std,
                        alphap=alphap_std,
                        i=i_best_schisto,Mp=Mp_mean,
                        xlab=c("E (activation energy)"),
                        plot_or_slope="plot",
                        color_or_black="color"){
  flux_E <- data.frame(E=E_r,flux=NA,slope=NA)
  j=2
  
  flux_E$flux[1] <- flux(Temp = Temp_std,E=E_r[1],
                         alphap=alphap,i=i,Mp=Mp)
  
  for(each_E in E_r[2:length(E_r)]){
    flux_E$flux[j] <- flux(E=E_r[j],Temp=Temp,
                           alphap=alphap,i=i,
                           Mp=Mp)
    flux_E$slope[j] <- (flux_E$flux[j-1] - flux_E$flux[j])/(flux_E$E[j-1]-flux_E$E[j])
    j <- j+1
  }

  flux_E_comp <- flux(E=base_E,Temp=Temp,alphap=alphap,i=i,Mp=Mp)
  flux_E$flux_diff <- (flux_E$flux-flux_E_comp)/(flux_E_comp)*100
  
  if(color_or_black=="color"){
    plot_flux_E_diff <- ggplot(flux_E,aes(x=E_r,
                                        y=flux_diff,
                                        color=abs(slope)))+
    theme_classic()+
    coord_cartesian(expand =FALSE)+
    geom_line()+
    geom_hline(yintercept = 0,color='gray')+ 
    ylim(c(-150,200))+
    theme(axis.title.y = element_blank(),
          legend.position = c(0.8,0.8)) + 
    xlab(xlab)+
    sc
    #scale_colour_gradientn(colours = myPalette(100), 
                           #limits=c(5e-13, 1e-11),
                           #oob=scales::squish)
  } else if(color_or_black=="black") {
    plot_flux_E_diff <- ggplot(flux_E,aes(x=E_r,
                                          y=flux_diff))+
      theme_classic()+
      coord_cartesian(expand =FALSE)+
      geom_line()+
      geom_hline(yintercept = 0,color='gray')+ 
      ylim(c(-150,200))+
      theme(axis.title.y = element_blank(),
            legend.position = c(0.8,0.8)) + 
      xlab(xlab)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
  
  df <- data.frame(min=min(flux_E$slope,na.rm=TRUE),
                   max=max(flux_E$slope,na.rm=TRUE))
  rownames(df) <- "E"
  
  # Output
  if (plot_or_slope == "plot"){
    plot_flux_E_diff
  } else if (plot_or_slope == "slope") {
    print(df)
  } else {
    print("Invalid entry - only accepts 'plot' or 'slope'")
  }
    }
plot_flux_E(plot_or_slope = "plot")
plot_flux_E(plot_or_slope = "plot",color_or_black = "black")

# Create jpg of multipanel figure
jpeg(file="analyses/sensitivity_analysis/flux_sensitivity_new.jpg",
     width = 10, height=2,units="in",res=300) # open jpg image
        plot_flux_alphap(color_or_black = "black")+theme(legend.position="none")+plot_flux_E(color_or_black = "black")+theme(legend.position="none")+ plot_flux_i(color_or_black = "black")+theme(legend.position="none")+plot_flux_Mp(color_or_black = "black")+theme(legend.position="none")+plot_flux_temp(color_or_black = "black")+theme(legend.position="none")+ plot_layout(ncol=5)+theme(legend.position="none")+plot_annotation(tag_levels = "A")

dev.off() # close jpg image

pdf(file="analyses/sensitivity_analysis/flux_sensitivity_new.pdf",
     width = 10, height=2) # open jpg image
plot_flux_alphap(color_or_black = "black")+theme(legend.position="none")+plot_flux_E(color_or_black = "black")+theme(legend.position="none")+ plot_flux_i(color_or_black = "black")+theme(legend.position="none")+plot_flux_Mp(color_or_black = "black")+theme(legend.position="none")+plot_flux_temp(color_or_black = "black")+theme(legend.position="none")+ plot_layout(ncol=5)+theme(legend.position="none")+plot_annotation(tag_levels = "A")
dev.off() # close image


jpeg(file="analyses/sensitivity_analysis/flux_sensitivity_2.jpg",width = 6, height=5.75,units="in",res=300) # open jpg image
plot_flux_alphap(Mp=Mp_min,xlab="",color_or_black = "black")+
  plot_flux_Mp(alphap=alphap_min,xlab="",color_or_black = "black")+
  plot_flux_temp(E=E_min,xlab="",color_or_black = "black")+
  plot_flux_alphap(Mp=Mp_mean,xlab="",color_or_black = "black")+
  plot_flux_Mp(alphap=alphap_std,xlab="",color_or_black = "black")+
  plot_flux_temp(E=E_std,xlab="",color_or_black = "black")+
  plot_flux_alphap(Mp=Mp_max,color_or_black = "black")+
  plot_flux_Mp(alphap=alphap_max,color_or_black = "black")+
  plot_flux_temp(E=E_max,color_or_black = "black")+
  plot_layout(guides="collect")+plot_annotation(tag_levels = "A")
dev.off() # close jpg image

