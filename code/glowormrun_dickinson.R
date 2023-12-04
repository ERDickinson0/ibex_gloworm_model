###############################################################################
#                         Gloworm simulations                                 #
#                                                                             #
#                       Author: Eleanor Dickinson                             #
#                       Date: November 2023                                   #
#  Code to run simulations of parasite transmission dynamics and              #
#  summarise results                                                          #
#  Simulations can be run with different host parameters                      #
#  Model outputs are saved as .csv files                                      #
###############################################################################

rm(list=ls())
setwd("~/PhD/Parasites/Glo worm")
source("glowormfunctions_dickinson.R")

###############################################################################
library(gridExtra)
library(ncdf4) # accessing the nc files from eobs climate data
library(dplyr)
###############################################################################

###############################################################################
#  Set up initial conditions                                                  #
###############################################################################
### Specify the location and time frame of the simulation ####
latitude <-  45.23 
longitude <-  6.07
ss <-  "2018-01-01" # start date
ee <-  "2018-12-31" # end date
compA <- 1 #specify the compartment for simulation
compB <- 2 #specify the compartment for simulation
compC <- 3 #specify the compartment for simulation
biomass <- (242.63 *10)/100*40 # Schweiger paper (ibex core foraging areas): 
#convert grams per m2 to kilograms per hectare and correct for dry herbage

### Loading ibex data ####
df <- read.csv("belledonnedata_FULL.csv") # parameters
df$date <- as.Date(df$date, "%d/%m/%Y") # ensure date is correct format
ibexage <- read.csv("levionazages.csv") # age data for ibex

### Calculate total density
df$density <- df$ibexdensity + df$sheepdensity
#-------------------------------------------------------------------------------

###############################################################################
#  Host calculations                                                          #
###############################################################################

##  Calculate FEC's when combining ibex and livestock - weighting the average 
### for the  density of each species
for(i in  1:nrow(df)){
  if(df$ibexdensity[i] == 0 & df$sheepdensity[i] == 0){
    df$totalweight[i] <- 0
  }else{
    df$propibex[i] <- ((df$ibexdensity[i]/
                          (df$ibexdensity[i]+df$sheepdensity[i]))*100)
    df$proplive[i] <- ((df$sheepdensity[i]/
                          (df$ibexdensity[i]+df$sheepdensity[i]))*100)
    df$FEC[i] <- (((df$FEC_ibexnew[i]*df$propibex[i])+
                     (df$FEC_sheep[i]*df$proplive[i]))/100)
  }
}

## Calculate weight when combining ibex and livestock - weighting the average 
# for the density of each species
for(i in  1:nrow(df)){
  if(df$ibexdensity[i] == 0 & df$sheepdensity[i] == 0){
    df$totalweight[i] <- 0
  }else{
    df$propibex[i] <- ((df$ibexdensity[i]/
                          (df$ibexdensity[i]+df$sheepdensity[i]))*100)
    df$proplive[i] <- ((df$sheepdensity[i]/
                          (df$ibexdensity[i]+df$sheepdensity[i]))*100)
    df$totalweight[i] <- (((df$weight_ibex[i]*df$propibex[i])+
                             (df$weight_sheep[i]*df$proplive[i]))/100)
  }
}

#-------------------------------------------------------------------------------

###############################################################################
#  Climate data                                                               #
#  We used the EOBS gridded dataset for Europe at 0.1 degree resolution       #
#  Daily timestep                                                             #
#  https://www.ecad.eu/download/ensembles/download.php                        #
#  (Other climate data could be incorporated)                                 #
###############################################################################

### Load the climate data using the E-OBS function
temp <-  eobs(data = "tg_ens_mean_0.25deg_reg_v20.0e.nc", var = "tg", 
              start = ss, end = ee, lat = latitude, lon = longitude)
precip <-  eobs(data = "rr_ens_mean_0.25deg_reg_v20.0e.nc", var = "rr", 
                start = ss, end = ee, lat = latitude, lon = longitude)
clim <- data.frame(temp,precip) # put both in a single data frame

### Calculate for each compartment - correct the temp and rain for the altitude
climA <- alt_correction(data = clim, compartment = compA)
climB <- alt_correction(data = clim, compartment = compB)
climC <- alt_correction(data = clim, compartment = compC)

### Plot the weather conditions 
climA$days <- seq(1, nrow(climA), by=1)
climB$days <- seq(1, nrow(climB), by=1)
climC$days <- seq(1, nrow(climC), by=1)
climate <- rbind(climA, climB, climC)

### Add season
climate$season <- NA
for(i in 1:nrow(climate)){
  if(climate$days[i] <= 59){
    climate$season[i] <- "Winter"
  }else if(climate$days[i] >= 335){
    climate$season[i] <- "Winter"
  }else if(climate$days[i] >= 60 && climate$days[i] <= 151){
    climate$season[i] <- "Spring"
  }else if(climate$days[i] >= 152 && climate$days[i] <= 243){
    climate$season[i] <- "Summer"
  }else{
    climate$season[i] <- "Autumn"
  }
}
#-------------------------------------------------------------------------------

###############################################################################
#  Running the simulation                                                     #
#                                                                             #
# - Density and FECS  arguments can be changed to reflect the target host     #
#   species                                                                   #
###############################################################################
#### COMPARTMENT A ####
### load initial values
initial_values <-  init.vals()

### Create an empty data frame
run <- data.frame(compartment = as.numeric(), 
                  rep = as.numeric(),
                  time = as.numeric(), 
                  E = as.numeric(), 
                  L = as.numeric(), 
                  L3f = as.numeric(), 
                  L3p = as.numeric(), 
                  Pasture = as.numeric(), 
                  Soil = as.numeric(), 
                  Herbage = as.numeric(), 
                  Herbperkg = as.numeric())

### Select the data for compartment A
ibexA <- subset(df, compartment == "A") 

### Run the simulation for 10 years 
runA <- run
w = 1
while(w<10){
  runA1 <- gloworm_full(init = initial_values, 
                        ss = ss, ee = ee, lat = latitude, 
                        temp = climA$temp2, precip = climA$precip2,
                        species = "teladorsagia", 
                        weight = ibexA$totalweight, ages = ibexage$x, 
                        ibexdates = ibexA$date, kgDM = biomass, 
                        density = ibexA$density,   fecs = ibexA$FEC)
  runA1$rep <- w
  runA <- rbind(runA, runA1)
  
  initial_values$L3faeces <- tail(runA$L3f, 1)
  initial_values$L3herbage <- tail(runA$Herbage, 1)
  w <- w+1
}

#### COMPARTMENT B ####
### load initial values
initial_values <-  init.vals()

### Select the data for compartment B
ibexB <- subset(df, compartment == "B") 

### Run the simulation for 10 years 
runB <- run
w = 1
while(w<10){
  runB1 <- gloworm_full(init = initial_values, 
                        ss = ss, ee = ee, lat = latitude, 
                        temp = climB$temp2, precip = climB$precip2,
                        species = "teladorsagia", 
                        weight = ibexB$totalweight, ages = ibexage$x, 
                        ibexdates = ibexB$date, kgDM = biomass, 
                        density = ibexB$density, fecs = ibexB$FEC)
  runB1$rep <- w
  runB <- rbind(runB, runB1)
  
  initial_values$L3faeces <- tail(runB$L3f, 1)
  initial_values$L3herbage <- tail(runB$Herbage, 1)
  w <- w+1
}

#### COMPARTMENT C ####
initial_values <-  init.vals()

### Select the data for compartment C
ibexC <- subset(df, compartment == "C") 

### Run the simulation for 10 years 
w = 1
runC <- run
while(w<10){
  runC1 <- gloworm_full(init = initial_values, 
                        ss = ss, ee = ee, lat = latitude, 
                        temp = climC$temp2, precip = climC$precip2, 
                        species = "teladorsagia", 
                        weight = ibexC$totalweight, ages = ibexage$x, 
                        ibexdates = ibexC$date, kgDM = biomass,
                        density = ibexC$density, fecs = ibexC$FEC)
  runC1$rep <- w
  runC <- rbind(runC, runC1)
  
  initial_values$L3faeces <- tail(runC$L3f, 1)
  initial_values$L3herbage <- tail(runC$Herbage, 1)
  w <- w+1
}

###############################################################################
# Compile and export the model output                                         #
###############################################################################

### compile the outputs
runA$compartment <- "A"
runB$compartment <- "B"
runC$compartment <- "C"
runA$Species <- "Alpine ibex and sheep"
runB$Species <- "Alpine ibex and sheep"
runC$Species <- "Alpine ibex and sheep"
fullrun <- rbind(runA, runB, runC)

### Export the output
write.csv(fullrun, file = "Modeloutput.csv")
#------------------------------------------------------------------------------