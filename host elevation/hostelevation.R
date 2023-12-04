###############################################################################
#                            Host elevation                                   #
#                                                                             #
#                       Author: Eleanor Dickinson                             #
#                       Date: November 2023                                   #
#                                                                             #
# 1. Host movement and elevation.                                             #
#                                                                             #
# 2. Predictors of ibex movement.                                             #
#                                                                             #
# 3. Predictors of sheep movement.                                            #
#                                                                             #
###############################################################################

rm(list=ls())
setwd("~/PhD/Parasites/BelledonneMovement")
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(MuMIn)
source("Altitude_function.R")
library(ncdf4) 
library(merTools)

data <- read.csv("avg_elev_week_2017-19.csv")
data$year <- as.factor(as.character(data$year))

###############################################################################
# 1. Host movement and elevation.                                             #
###############################################################################

#### Ibex elevation ####
ibex <- subset(data, species == "ibex")
mod <- lm(elevation ~ year + season, data = ibex)
summary(mod)

model <- lmer(elevation ~ year + (1|id:season), data = ibex)

#### Sheep presence ####
sheep <- subset(data, species == "sheep")
model3 <- lmer(elevation ~ year + week + (1|id), data = sheep)


### label compartments ###
newdata <- data 
for(i in 1:nrow(newdata)) {
  if(newdata$elevation[i] < 1800){
    newdata$compartment[i] <- "A"
  }else if(newdata$elevation[i] > 2150){
    newdata$compartment[i] <- "C"
  }else{
    newdata$compartment[i] <- "B"
  }
}

data2 <- newdata %>% 
  group_by(compartment) %>% 
  summarise (mean = mean(elevation),median = median(elevation), 
             sd = sd(elevation))

#### Calculate the number of observations in each compartment ####
newdata$obsv <- 1

df1 <- ddply(newdata, c("week", "species", "compartment"),  
                 function(x) 
                 {
                   cv <- sum(x$obsv)
                   
                   data.frame(num = cv)
                 })
df1

df2 <- ddply(newdata, c("week", "species"),  
                  function(x) 
                  {
                    cv <- sum(x$obsv)
                    
                    data.frame(total = cv)
                  })
df2
####### Calculate the proportion of the population observed in each compartment, 
####### each week

weeks <- seq(1, 53, 1)

#### Compartment A, ibex 
Ibex.A <- data.frame(Week = weeks)

Ibex.A$Species <- "Ibex"
Ibex.A$Compartment <- "A"

Ibex.A$No <- as.numeric(paste(df1[df1$species == "ibex" & 
                               df1$compartment == "A",]$num[match(Ibex.A$Week,
                  df1[df1$species == "ibex" & df1$compartment == "A",]$week)]))

Ibex.A[is.na(Ibex.A$No),]$No <- 0

Ibex.A$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.A$Week,
                                  df2[df2$species == "ibex",]$week)]))
head(Ibex.A)

#### Compartment B, ibex 
Ibex.B <- data.frame(Week = weeks)

Ibex.B$Species <- "Ibex"
Ibex.B$Compartment <- "B"

Ibex.B$No <- as.numeric(paste(df1[df1$species == "ibex" & df1$compartment == "B",]$num[match(Ibex.B$Week,
                     df1[df1$species == "ibex" & df1$compartment == "B",]$week)]))

Ibex.B[is.na(Ibex.B$No),]$No <- 0

Ibex.B$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.B$Week,
                               df2[df2$species == "ibex",]$week)]))
head(Ibex.B)

#### Compartment C, ibex 
Ibex.C <- data.frame(Week = weeks)

Ibex.C$Species <- "Ibex"
Ibex.C$Compartment <- "C"

Ibex.C$No <- as.numeric(paste(df1[df1$species == "ibex" & df1$compartment == "C",]$num[match(Ibex.B$Week,
                 df1[df1$species == "ibex" & df1$compartment == "C",]$week)]))

Ibex.C[is.na(Ibex.C$No),]$No <- 0

Ibex.C$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.C$Week,
                              df2[df2$species == "ibex",]$week)]))
head(Ibex.C)

#### Compartment A, sheep
sheep.A <- data.frame(Week = weeks)

sheep.A$Species <- "sheep"
sheep.A$Compartment <- "A"

sheep.A$No <- as.numeric(paste(df1[df1$species == "sheep" & df1$compartment == "A",]$num[match(sheep.A$Week,
                 df1[df1$species == "sheep" & df1$compartment == "A",]$week)]))

sheep.A[is.na(sheep.A$No),]$No <- 0

sheep.A$Tot <- as.numeric(paste(df2[df2$species == "sheep",]$total[match(sheep.A$Week,
                         df2[df2$species == "sheep",]$week)]))
sheep.A[is.na(sheep.A$Tot),]$Tot <- 0

head(sheep.A)

#### Compartment , sheep
sheep.B <- data.frame(Week = weeks)

sheep.B$Species <- "sheep"
sheep.B$Compartment <- "B"

sheep.B$No <- as.numeric(paste(df1[df1$species == "sheep" & df1$compartment == "B",]$num[match(sheep.B$Week,
                  df1[df1$species == "sheep" & df1$compartment == "B",]$week)]))

sheep.B[is.na(sheep.B$No),]$No <- 0

sheep.B$Tot <- as.numeric(paste(df2[df2$species == "sheep",]$total[match(sheep.B$Week,
                                       df2[df2$species == "sheep",]$week)]))
sheep.B[is.na(sheep.B$Tot),]$Tot <- 0

head(sheep.B)

#### Compartment C, sheep
sheep.C <- data.frame(Week = weeks)

sheep.C$Species <- "sheep"
sheep.C$Compartment <- "C"

sheep.C$No <- as.numeric(paste(df1[df1$species == "sheep" & df1$compartment == "C",]$num[match(sheep.B$Week,
                 df1[df1$species == "sheep" & df1$compartment == "C",]$week)]))

sheep.C[is.na(sheep.C$No),]$No <- 0

sheep.C$Tot <- as.numeric(paste(df2[df2$species == "sheep",]$total[match(sheep.C$Week,
                                         df2[df2$species == "sheep",]$week)]))
sheep.C[is.na(sheep.C$Tot),]$Tot <- 0

head(sheep.C)


#### Combine and calculate proportion
newdf <- rbind(Ibex.A, Ibex.B, Ibex.C, sheep.A, sheep.B, sheep.C)

newdf2 <- ddply(newdf, c("Week", "Species", "Compartment"),  
             function(x) 
             {
               n <- x$No
               t <- x$Tot
               cv <- x$No/x$Tot
               
               data.frame(No = n, Tot = t, prop = cv)
             })

newdf2[is.na(newdf2$prop),]$prop <- 0
head(newdf2)

write.csv(newdf2, file = "belledonne_proportionsALL.csv")

##### Add to the dataframe #############

#-------------------------------------------------------------------------------

###############################################################################
# 3. Predictors of host movement                                              #
###############################################################################

###### Climate conditions #####
latitude <-  45.23 # location of simulation 
longitude <-  6.07
ss <-  "2017-05-21" # 2017 week 20
ee <-  "2019-05-27" # 2019 week 15

#### load the data 
temp <-  eobs(data = "tg_ens_mean_0.25deg_reg_v20.0e.nc", var = "tg", 
              start = ss, end = ee, lat = latitude, lon = longitude)
precip <-  eobs(data = "rr_ens_mean_0.25deg_reg_v20.0e.nc", var = "rr", 
                start = ss, end = ee, lat = latitude, lon = longitude)
clim <- data.frame(temp,precip) # put both in a single data frame

# Mean weekly values
date <- seq(from = as.Date(ss), to = as.Date(ee), by = 1)
clim$date <- date
clim$week <- strftime(clim$date, format = "%V")
clim$year <- strftime(clim$date, format = "%Y")

### average weekly temperatures
clim2 <- ddply(clim, c("year", "week"),  
               function(x) 
               {
                 ab <- mean(x$temp)
                 bb <- sum(x$precip)
                 
                 data.frame(temp = ab, precip = bb)
               })
clim2$week <- as.numeric(clim2$week)
clim2$key <- paste(clim2$year, clim2$week, sep = "-")

###
ibex$week <- as.numeric(ibex$week)
ibex$key <- paste(ibex$year, ibex$week, sep = "-")
ibex$year <- as.factor(ibex$year)
head(ibex)

#### Add climate to host movement and altitude correct ####
ibex$temp <- as.numeric(paste(clim2$temp[match(ibex$key, clim2$key)]))
ibex$precip <- as.numeric(paste(clim2$precip[match(ibex$key, clim2$key)]))
ibex <- na.omit(ibex)

ibex$temp2 <-  as.numeric(ibex$temp - ((ibex$elevation - 1839)*(0.0047)))
for(i in 1:nrow(ibex)){
  if(ibex$precip[i] <= 0){
    ibex$precip2[i] <- 0
  }else{
    ibex$precip2[i] <- as.numeric(ibex$precip[i] + ((ibex$elevation[i] - 1839)*
                                                      (3.2 *10^(-4))))
  }}


##### Effect of temperature etc..
model <- lmer(elevation ~ year + temp  + precip + season + (1|id:site),
              data = ibex)

newX <- expand.grid(temp = seq(from = -3, to = 23, length = 27),
                    season = levels(ibex$season),
                    id= levels(ibex$id))
newY <- predict(model, newdata = newX, re.form = NA)
fm1W <- confint(model, method="Wald")
head(newY)
newX$elevation <- newY
newX$no <- rep(1:27, each = 1, times = 360, len = 2808)
newX$no <- as.factor(newX$no)

pred <- cbind(newX, predictInterval(model1, newX))

predicted <- pred %>% 
  dplyr:::group_by(no) %>% 
  dplyr:::summarise(temp = mean(temp), elevation = mean(fit), 
                    upr = mean(upr), lwr = mean(lwr))

#### Avergae over three years
ibex2 <- ibex %>% 
  dplyr:::group_by(week) %>% 
  dplyr:::summarise(temp = mean(temp), precip = mean(precip), temp2 = mean(temp2), 
                    precip2 = mean(precip2), elevation = mean(elevation))
head(ibex2)

#### Ibex will be 24.67m higher per 1 degree C temperature increase ####
#------------------------------------------------------------------------------

###############################################################################
#### Projected temperature 
# Elevation of ibex each week
# temperature each week of elevation
# mean projected temperature per week 
# calculate projected mean weekly elevation 
# use that to calculate the proportion of ibex in each compartment
###############################################################################


#### Predict elevation on projected climate data
proj <- read.csv("RCP85_2065-2095_365.csv")
head(proj)
# Mean weekly values
str(proj)
proj$Date <- as.Date(proj$Date, format = "%d/%m/%Y")
proj$week <- strftime(proj$Date, format = "%V")
proj$year <- strftime(proj$Date, format = "%Y")

### average weekly temperatures
proj2 <- ddply(proj, c("week"),  
               function(x) 
               {
                 ab <- mean(x$Temperature)
                 bb <- sum(x$Preciptiation)
                 
                 data.frame(temp = ab, precip = bb)
               })
proj2$week <- as.numeric(proj2$week)

### add proj temp to ibex data
ibex2$projtemp <- as.numeric(paste(proj2$temp[match(ibex2$week,proj2$week)]))
ibex2$projprecip <- as.numeric(paste(proj2$precip[match(ibex2$week,proj2$week)]))

### adjust for elevation
ibex2$projtemp2 <-  as.numeric(ibex2$projtemp - ((ibex2$elevation - 1839)*
                                                   (0.0047)))
ibex2$projprecip2 <- as.numeric(0)
for(i in 1:nrow(ibex2)){
  if(ibex2$projprecip[i] <= 0){
    ibex2$projprecip2[i] <- 0
  }else{
    ibex2$projprecip2[i] <- as.numeric(ibex2$projprecip[i] + 
                                         ((ibex2$elevation[i] - 1839)*
                                            (3.2 *10^(-4))))
  }}

### project ibex elevation with climate change 
# put the temperature back on the full data 
ibex$tempproj <- as.numeric(paste(ibex2$projtemp2[match(ibex$week,ibex2$week)]))

#### project elevation
x <- ibex$tempproj - ibex$temp2
y <- ibex$elevation + (24.67 * x)
ibex$elevationproj <- y

############## Projected proportion of ibex in each compartment ############
##### label compartments ###
newdata <- ibex
newdata$compartment <- NA

for(i in 1:nrow(newdata)) {
  if(newdata$elevationproj[i] < 1800){
    newdata$compartment[i] <- "A"
  }else if(newdata$elevationproj[i] > 2150){
    newdata$compartment[i] <- "C"
  }else{
    newdata$compartment[i] <- "B"
  }
}

newdata$obsv <- 1

df1 <- ddply(newdata, c("week", "species", "compartment"),  
             function(x) 
             {
               cv <- sum(x$obsv)
               
               data.frame(num = cv)
             })
head(df1)

df2 <- ddply(newdata, c("week", "species"),  
             function(x) 
             {
               cv <- sum(x$obsv)
               
               data.frame(total = cv)
             })
head(df2)

####### Calculating porportions
weeks <- seq(1, 52, 1)

#### Compartment A
Ibex.A <- data.frame(Week = weeks)

Ibex.A$Species <- "Ibex"
Ibex.A$Compartment <- "A"

Ibex.A$No <- as.numeric(paste(df1[df1$species == "ibex" & df1$compartment == "A",]$num[match(Ibex.A$Week,
                        df1[df1$species == "ibex" & df1$compartment == "A",]$week)]))

Ibex.A[is.na(Ibex.A$No),]$No <- 0

Ibex.A$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.A$Week,
                          df2[df2$species == "ibex",]$week)]))
head(Ibex.A)

#### Compartment B
Ibex.B <- data.frame(Week = weeks)

Ibex.B$Species <- "Ibex"
Ibex.B$Compartment <- "B"

Ibex.B$No <- as.numeric(paste(df1[df1$species == "ibex" & df1$compartment == "B",]$num[match(Ibex.B$Week,
                        df1[df1$species == "ibex" & df1$compartment == "B",]$week)]))

Ibex.B[is.na(Ibex.B$No),]$No <- 0

Ibex.B$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.B$Week,
                         df2[df2$species == "ibex",]$week)]))
head(Ibex.B)

#### Compartment C
Ibex.C <- data.frame(Week = weeks)

Ibex.C$Species <- "Ibex"
Ibex.C$Compartment <- "C"

Ibex.C$No <- as.numeric(paste(df1[df1$species == "ibex" & df1$compartment == "C",]$num[match(Ibex.B$Week,
                        df1[df1$species == "ibex" & df1$compartment == "C",]$week)]))

Ibex.C[is.na(Ibex.C$No),]$No <- 0

Ibex.C$Tot <- as.numeric(paste(df2[df2$species == "ibex",]$total[match(Ibex.C$Week,
                         df2[df2$species == "ibex",]$week)]))
head(Ibex.C)


#### Combine and calculate proportion
newdf <- rbind(Ibex.A, Ibex.B, Ibex.C)

str(newdf)

newdf2 <- ddply(newdf, c("Week", "Species", "Compartment"),  
                function(x) 
                {
                  n <- x$No
                  t <- x$Tot
                  cv <- x$No/x$Tot
                  
                  data.frame(No = n, Tot = t, prop = cv)
                })

newdf2[is.na(newdf2$prop),]$prop <- 0
head(newdf2)
#------------------------------------------------------------------------------

###############################################################################
# Predictors of sheep movement                                                #
###############################################################################

#### Host movement data ####
data <- read.csv("avg_elev_week_alldata.csv")

head(data)
str(data)

data$week <- as.numeric(data$week)
data$key <- paste(data$year, data$week, sep = "-")
data$year <- as.factor(data$year)
data$site <- as.factor(data$site)
sheep <- subset(data, species == "sheep")
clim2$elevation <- as.numeric(paste(sheep$elevation[match(clim2$key, 
                                                          sheep$key)]))
clim2$site <- paste(sheep$site[match(clim2$key, sheep$key)])

### Add climate data 
#### sheep df to model elevation, clim2 df to model presence/absence
clim2$sheep <- ifelse(is.na(clim2$elevation), "0", "1")
clim2$sheep <- as.numeric(clim2$sheep)
clim2$year <- as.numeric(clim2$year)

#### Model sheep presence/absence
#### Predict the probabilities
newX <- expand.grid(Temp = seq(from = min(clim2$Temp), to = max(clim2$Temp), length = 15))

# predict(model, new.df) will just give you the logit (z) values
# We want the probabilities so use the "response" type
probs <- predict(model3, newX, "response")
newdata <- cbind(newX, probs)
head(newdata)
plot(clim2$Temp, clim2$sheep, pch = 16, xlab = "Temp", ylab = "Sheep presence")
lines(newdata$Temp, newdata$probs)

#### Predict elevation on projected climate data
proj <- read.csv("RCP85_2065-2095_365.csv")
head(proj)
# Mean weekly values
str(proj)
proj$Date <- as.Date(proj$Date, format = "%d/%m/%Y")
proj$week <- strftime(proj$Date, format = "%V")
proj$year <- strftime(proj$Date, format = "%Y")

### average weekly temperatures
proj2 <- proj %>% 
  dplyr:::group_by(year, week) %>% 
  dplyr:::summarize(Temp = mean(Temperature), Precip = mean(Preciptiation))

proj2$week <- as.numeric(proj2$week)
proj2$key <- paste(proj2$year, proj2$week, sep = "-")
head(proj2)

proj2$sheep <- ifelse(proj2$Temp >= 15.55482, "P", "A")
#------------------------------------------------------------------------------
