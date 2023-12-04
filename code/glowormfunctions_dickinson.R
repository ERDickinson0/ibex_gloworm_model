###############################################################################
#                         Gloworm simulations                                 #
#                                                                             #
#                       Author: Eleanor Dickinson                             #
#                       Date: November 2023                                   #
#                                                                             #
#### This script contains the following functions for the simulations:        #
#### 1. Parasite Date interpolation function                                  #
#### 2. Climate data function to extract the climate data (EOBS)              #
#### 3. Altitude correction function for climate data                         #
#### 4. Initial values function                                               #
#### 5. Dry matter intake function                                            #
#### 6. Faeces production function                                            #
#### 7. Glowowrm function                                                     #
###############################################################################

setwd("~/PhD/Parasites/Glo worm")

### Load up the relevant packages. ###
library(geosphere)
library(forecast)
library(deSolve)
library(chron)
library(plotrix)
library(imputeTS)
library(forecast)
library(scales)
library(tidyverse)

################# Parasite Date interpolation function ########################
# Combines the egg input and the timeframe you have specified 
#### Arguments: - para = list of numbers of faecal egg counts
####            - dates = list of dates for each FEC including start date and 
####                      end date.
#### Ouput: A list of parasite counts for each day within the specified period

parainterp <-  function(para, dates, method="linear") {
  indices <-  which(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days") %in% as.Date(dates))
  days <-  seq(1, length(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days")), 1)
  interpfun <-  approxfun(y = para, x = indices, method = method)
  
  return(interpfun(days))
}


############# E-OBS function #################################################
#### Climate data function to extract the climate data from the eobs files ###
#### Arguments: - data = the data file with the climatic data 
####          - var = the variable you want to extract eg. temp
####          - lat = latitude of the study site
####          - lon = longitude of the study site
####          - start = beginning of the period for simulation
####          - end = end of the simulation period
#### Output = a string of numbers denoting the variable of interested daily for 
####         the simulation period

eobs <-  function(data, var, lat, lon, start, end) {
  
  dat <-  nc_open(data)
  
  latitude <- ncvar_get(dat,"latitude")
  longitude <- ncvar_get(dat,"longitude")  
  
  t <- ncvar_get(dat,"time")
  tunits <- ncatt_get(dat,"time","units")
  
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <-  as.integer(unlist(tdstr)[2])
  tday <-  as.integer(unlist(tdstr)[3])
  tyear <-  as.integer(unlist(tdstr)[1])
  tchron <-  as.Date(chron(t, origin = c(tmonth, tday, tyear)))
  
  start. <-  c(which.min(longitude < lon)-1,which.min(latitude < lat)-1,which(tchron == start))
  count <-  c(1,1,which(tchron == end)+1-which(tchron == start))
  
  dat2 <- ncvar_get(dat, var, start., count)
  
  nc_close(dat)
  
  return(dat2)
}


#### Altitude correction function ####
### This function corrects the temperature and rainfall datafor altitude
### It is specific for Alp compartments labelled A, B, C
### Numbers in fuction (first model) = 1682 = altitude of EOBS dataset 
###                     2000 (A), 2350 (B), 2700 (C) mean altitude of compartments  
### Arguments: - data = a dataframe with two columns "temp" and "precip"
alt_correction <- function(data, compartment){
  for(i in 1:nrow(data)){
    if (compartment == 1){
      data$temp2[i] <-  as.numeric(data$temp[i] - ((1661 - 1839)*(0.0047)))
    }else if(compartment == 2){
      data$temp2[i]  <- as.numeric(data$temp[i] - ((1915 - 1839)*(0.0047)))
    }else if(compartment == 3){
      data$temp2[i] <- as.numeric(data$temp[i] - ((2300 - 1839)*(0.0047)))
    }else{
      NA
    }
  }
  
  for(i in 1:nrow(data)){
    if(data$precip[i] == 0){
      data$precip2[i] <- 0
    }else{
      if (compartment == 1){
        data$precip2[i] <- as.numeric(data$precip[i] + ((1661 - 1839)*(3.2 *10^(-4))))
      }else if(compartment == 2){
        data$precip2[i] <- as.numeric(data$precip[i] + ((1915 - 1839)*(3.2 *10^(-4))))
      }else if(compartment == 3){
        data$precip2[i] <- as.numeric(data$precip[i] + ((2300 - 1839)*(3.2 *10^(-4))))
      }else{
        NA
      }
    }
  }
  data$compartment <- compartment
  return(data)
}


########################## Initial values function ###########################
# Function to generate default initial values - the situation at the #
# start of the season - change if you know how many larvae are on the #
# pasture at the start of the season #

init.vals <-  function(Eggs =0, L1L2 =0, L3faeces=0, L3herbage=0) {
  
  v.mig <-  (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * (temp^2))))
  L3herbage <-  L3herbage * (1/v.mig[1])
  
  .init <-  cbind(Eggs, L1L2, L3faeces, L3herbage)
  return(as.data.frame(.init))
}

###################### Dry matter intake function ############################
# Cattle and sheep dry matter intake calculator - I assume ibex are the same as sheep
# Arguments: age = age of group
#            weight = live weight of group
#            host species = host species (options atm cow/sheep)

dmi <-  function(lwt, age, host_species, milk = 0) {
  
  if(host_species==1) { #sheep (or ibex)
    DMI = 0.0545*(lwt^0.86)
    if(any(age<51)) {
      interp_DMI = approxfun(y = c(0, DMI[which(age==50)[1]]), x = c(1, 50), method = 'linear')
      DMI[which(age>0&age<=50)] = interp_DMI(age[which(age>0&age<=50)])
    }}
  
  if(host_species==2) { #cattle
    DMI = rep(NA, length(lwt))
    for(i in 1:length(lwt)) {
      if(lwt[i]<=400) DMI[i]=(31.4-0.03*lwt[i])*(lwt[i]/1000) else (DMI[i] = 0.025*lwt[i]+0.1*milk)
    }}
  
  return(DMI)
  
}

###################### Faeces production function ############################
# Cattle and Sheep faeces production calculator - what will this be for ibex? 
# * find out what the number in this function are for to know what i need for ibex
# need to change this to ibex and sheep.
# 1 = sheep, 2 = cattle

faeces <-  function(lwt, host_species) {
  if(host_species==1) { f = lwt*7}
  if(host_species==2) {
    f = ifelse(lwt<=250, (lwt*0.0811)*(1-0.34), (lwt*0.0181+17.8)*(1-0.34))
    f = f*1000
  }
  return(f)
}


##################### GLO WORM movement function #########################
#### This function
#### Arguments: t =
####            y = 
#### output: 

gloworm_full <-  function (init, ss, ee, lat, temp, precip, species, # parameters for parasites
                           weight, ages, ibexdates, kgDM, # parameters for dmi
                           dates, density, fecs) {
  
  date.range <- seq(as.Date(ss), as.Date(ee), "days") # set the time period of the study
  if (length(date.range) != length(temp)) # checks that the climate data and dates match 
    stop("Error: length of climatic data and dates do not match. Ensure climatic data 
         correspond to the start and end dates (one entry per day)")
  
  global.t <- seq(1, length(date.range))  # length of study
  photoperiod <- daylength(lat = lat, doy = date.range) # days lengths for time period
  
  # Specify the parasite species
  if (species == "teladorsagia") { ## THIS SECTION SPECIFIES THE PARAMETERS FOR TELADORSAGIA CIRCUMCINCTA ## 
    Pt = (4.95 * exp(0.062 * temp))/100
    PET = 0.55 * ((photoperiod/12)^2) * Pt * 25.4
    PET7 = NA
    PETm7 = NA
    Precip7 = NA
    Precipm7 = NA
    for (z in 1:length(PET)) {
      min = z
      max = ifelse((z + 7) > length(PET), length(PET), # 7 means 7 days of no rain and thats when it wont develop
                   z + 7)
      PET7[z] = sum(x = PET[min:max])
      Precip7[z] = sum(x = precip[min:max])
      max = z
      min = ifelse((z - 7) < 1, 1, z - 7)
      PETm7[z] = sum(x = PET[min:max])
      Precipm7[z] = sum(x = precip[min:max])
      P.E.1 = Precip7/PET7
      P.E.2 = Precipm7/PETm7
    }
    dev.1 = pmax(0, -0.02085 + 0.00467 * temp)
    dev.1 = pmin(1, dev.1)
    mu.1 = pmin(1, exp(-1.62026 - 0.17771 * temp + 0.00629 * 
                         (temp^2)))
    mu.2 = mu.1
    mu.4 = pmin(1, exp(-4.58817 - 0.13996 * temp + 0.00461 * 
                         (temp^2)))
    mu.3 = pmin(1, 10 * mu.4)
    mu.5 = mu.3
    h.mig = ifelse(precip >= 2, 0.21, ifelse(P.E.2 >= 1, 
                                             0.025, 0))
    v.mig = (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * 
                           (temp^2))))
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = ifelse(P.E.1 < 1, 0.1, 1)
    init[1] = init[1] * egg.correction[1]
    print("Teladorsagia circumcincta parameters loaded successfully")
    
  }
  if (species == "ostertagia") { ## THIS SECTION SPECIFIES THE PARAMETERS FOR OSTERTAGIA OSTERTAGI ## 
    dev.1 = pmax(0, -0.07258 + 0.00976 * temp)
    dev.1 = pmin(1, dev.1)
    mu.1 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * 
                         (temp^2)))
    mu.2 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * 
                         (temp^2)))
    mu.4 = pmin(1, exp(-6.388 - 0.2681 * temp + 0.01633 * 
                         (temp^2) - 0.00016 * (temp^3)))
    mu.3 = pmin(1, mu.4 * 10)
    mu.5 = mu.3
    h.mig = ifelse(precip >= 2, 0.06, 0)
    v.mig = (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * 
                           (temp^2))))
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    init[5] <-  (init[4] * (1/v.mig[1])) - init[4]
    egg.correction <-  rep(1, length(global.t))
    print("Ostertagia ostertagi parameters loaded successfully")
  }
  # check the right species name
  if (species != "teladorsagia" & species != "ostertagia") 
    stop("Error: 'species' parameter undefined. Please define nematode species or 
         check spelling") 

  ## Specify the DMI 
  age.vector <- (mean(ages)*365) # specify the average age of the group

  # DMI function
  DMI <- parainterp(para = dmi(lwt = weight, age = age.vector, host_species = 1),
                    dates = ibexdates, method = "linear") 
return(DMI)
  # Estimate the proportion of available herbage consumed per day, 
  #and the daily instantaneous rate of intake 
  pDMI <-  (DMI/kgDM)
  rDMI <-  -log(1-pDMI)
  # Create inerpolation function for dry matter intake rate
  dry_matter_intake <-  approxfun(rDMI, method = "linear")
  if (length(date.range) != length(DMI)) # checks that the climate data and dates match 
    stop("Error: in DMI calculation")
  print("DMI parameters loaded successfully")
  
  ## Specify the group density
  SR <-  parainterp(para = density, dates = ibexdates, method = "constant") #stockingrate
  density_group <-  approxfun(SR, method = "constant") #interpolate density
  print("Stocking rate parameters loaded successfully")
  
  ## Create FEC values for egg input into system
  fec <-  parainterp(para = fecs, dates = ibexdates) # interpolate FEC's
  av.faeces <- faeces(lwt = weight, host_species = 1) # average fecal output 
  faeces_ha <- density * av.faeces# faeces per hec = ibex density * faeces per indiv
  faeces_ha2 <-  parainterp(para = faeces_ha, dates = ibexdates) # ibex density * faeces 
  FECs_total <-  fec*faeces_ha2 ## final object for the model ##
  print("FEC parameters loaded successfully")
  
  ## Assign the movement function
  gloworm_movement <-  function (t, para.init, para.par) {
    with(as.list(c(para.init, para.par)), {
      
      dev1 <-  dev1rate(t)
      mu1 <-  mu1rate(t)
      mu2 <-  mu2rate(t)
      mu3 <-  mu3rate(t)
      mu4 <-  mu4rate(t)
      mu5 <-  mu5rate(t)
      m1 <-  hmigrate(t)
      m2 <-  vmigrate(t)
      
      i <-  dry_matter_intake(t)      
      stock.rate <-  density_group(t) 
      g <-  graze(t)               # grazing paddock for t = timeperiod
      
      #equations for the parasites for the pasture
      dE = -(dev1 * 2 + mu1) * E
      dL = -(dev1 * 2 + mu2) * L + (dev1 * 2) * E
      dL3f = -(mu3 + m1) * L3f + (dev1 * 2) * L
      dL3p = -mu4 * (L3p * (1 - m2)) - mu5 * (L3p * m2) + m1 * L3f - (L3p * m2)*
        i*stock.rate*g # pasture as a whole including dry matter intake and stocking rate
      
      return(list(c(dE = dE, dL = dL, dL3f = dL3f, dL3p = dL3p)))
    })
  }
  print("model function loaded successfully")
  
  # load initial values
  para.init <-  c(E = initial_values$Eggs, 
                  L = initial_values$L1L2, 
                  L3f = initial_values$L3faeces, 
                  L3p = initial_values$L3herbage * (1/v.mig[1]))
  print("initial conditions for state variables loaded successfully")
  
  # Specify grazing in area
  g <- rep(1, length(date.range))
  graze <-  approxfun(ifelse(g==1, 1, 0), method = "constant") 
  var_event <-  rep("E", length(g))
  
  ## Create event data frame
  eventdat <-  data.frame(var = var_event, # a list of movement events
                          time = seq(1, length(global.t)), # each day for length of simulation
                          value = FECs_total, method = "add")
  print("Event data parameters loaded successfully")
  
  para.sol <-  lsoda(y = para.init, # user defined initial values 
                times = global.t, 
                func = gloworm_movement, 
                parms = NULL,
                events = list(data = eventdat))
  print("simulation successfull - check for warnings")
  
  para.sol <- data.frame(para.sol)# make the output a dataframe
  Pasture <- para.sol[, "L3p"]
  Soil <-  Pasture * (1 - (vmigrate(global.t))) # calculate the L3 in soil
  Herbage <-  Pasture * vmigrate(global.t) # calculte the L3 on herbage per hectare
  Herbperkg <- Herbage/kgDM # calculate the L£ per kg herbage
  para.sol <-  cbind(para.sol, Pasture, Soil, Herbage, Herbperkg)
  print("post-hoc calculations successful")
 
  para.sol$date <- as.Date(date.range, "%Y-%m-%d")
  print("Output dataframe with dates and posthoc calculations")
  return(para.sol)
  
}