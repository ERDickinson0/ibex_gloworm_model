###################### Running gloworm more efficiently ########################

rm(list=ls())
setwd("~/PhD/Parasites/Glo worm")
source("glowormfunctions_dickinson.R")
library(gridExtra)
library(ncdf4) # accessing the nc files from eobs
library(dplyr)

### Specify the location and time frame of the simulation ####
latitude <-  45.23 # location of simulation 
longitude <-  6.07
ss <-  "2018-01-01" # start date of simulation
ee <-  "2018-12-31" # end date
compA <- 1 #specify the compartment for simulation
compB <- 2 #specify the compartment for simulation
compC <- 3 #specify the compartment for simulation
biomass <- (242.63 *10)/100*40 # Schweiger paper (ibex core foraging areas): convert:
#convert grams per m2 to kilograms per hectare and correct for dry herbage

##### Loading ibex data ####
df <- read.csv("belledonnedata_FULL.csv")
str(df)
ibexage <- read.csv("levionazages.csv")
df$date <- as.Date(df$date, "%d/%m/%Y") # ensure date is correct format

#### Calculate total density
df$density <- df$ibexdensity + df$sheepdensity
head(df)

hist(df$density)
hist(df$FEC)

## Calculate FEC's when combining ibex and livestock - weighting the average for the 
# density of each species
#df$FEC <- df$FEC_ibex + df$FEC_sheep

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

## Calculate weight when combining ibex and livestock - weighting the average for the 
# density of each species
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

str(df)
sum <- df %>% 
  group_by(compartment, season) %>% 
  summarise(ibex=(FEC_ibexnew), sheep=(FEC_sheep))
sum
write.csv(sum, file = "meanFECS.csv")

##### Load the climate data using the E-OBS function ####
# using eobs function want to open the doc, 
# select the variable, start date, end date, lat, long.
temp <-  eobs(data = "tg_ens_mean_0.25deg_reg_v20.0e.nc", var = "tg", 
              start = ss, end = ee, lat = latitude, lon = longitude)
precip <-  eobs(data = "rr_ens_mean_0.25deg_reg_v20.0e.nc", var = "rr", 
                start = ss, end = ee, lat = latitude, lon = longitude)
clim <- data.frame(temp,precip) # put both in a single data frame

hist(clim$temp)
hist(clim$precip)
median(clim$precip)

## Calculate for each compartment - correct the temp and rain for the altitude
climA <- alt_correction(data = clim, compartment = compA)
climB <- alt_correction(data = clim, compartment = compB)
climC <- alt_correction(data = clim, compartment = compC)

#### Plot the weather conditions ####
climA$days <- seq(1, nrow(climA), by=1)
climB$days <- seq(1, nrow(climB), by=1)
climC$days <- seq(1, nrow(climC), by=1)
climate <- rbind(climA, climB, climC)

head(climate)

sumclim <- climate %>% 
  group_by(compartment) %>% 
  summarize(meanT = mean(temp2), sdT = sd(temp2),minT = min(temp2), maxT = max(temp2),
            meanP = median(precip2), iqr = IQR(precip2))
sumclim


### Add season ###
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

sum <- climate %>% 
  dplyr:::group_by(season, compartment) %>% 
  dplyr:::summarise(mTemp = mean(temp2), MPrecip = mean(precip2))
write.csv(sum, file = "sum.csv")

# plot
tp <- ggplot(climate, aes(x= days, y = temp))+
  geom_line(linetype = "dashed", size=0.6)+
  geom_line(aes(y =temp2, colour = compartment), size=0.6)+
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+
  scale_colour_manual(values= c("#28346E", "#364696", "#5F70C5"))+
  ggtitle( "A) Temperature")+ylab("Temperature (oC)")+ xlab("Month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90),
                   legend.position = 'none' )

rn <- ggplot(climate, aes(x= days, y = log(precip+1)), size=0.6)+
  geom_line(linetype = "dashed")+
  geom_line(aes(y =precip2, colour = compartment), size=0.6)+
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+
  ggtitle( "B) Precipitation") +ylab("Precipitation (mm)")+ xlab("Month")+
  scale_colour_manual(values= c("#28346E", "#364696", "#5F70C5"))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90),
                   legend.position = 'none')

#png(file="Belledonne_TempPrecip.png",width=1000,height=500,res=110)
grid.arrange(tp, rn, ncol = 2)
#dev.off()


### load intial values
initial_values <-  init.vals()

ibexA <- subset(df, compartment == "A") # subset the data by compartment
str(ibexA)
head(climB)

## Run the simulation ##

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

## Run the simulation ##
w = 1
runA <- run

while(w<10){
  runA1 <- gloworm_full(init = initial_values, # initial values
             ss = ss, ee = ee, lat = latitude, # start and end dates, and location
             temp = climA$temp2, precip = climA$precip2, # temp and precip parameters
             species = "teladorsagia", # parameter for parasite species
             weight = ibexA$totalweight, ages = ibexage$x, 
             ibexdates = ibexA$date, kgDM = biomass, # parameters for dmi
             density = ibexA$density,   fecs = ibexA$FEC)
runA1$rep <- w
runA <- rbind(runA, runA1)

initial_values$L3faeces <- tail(runA$L3f, 1)
initial_values$L3herbage <- tail(runA$Herbage, 1)
w <- w+1
}

head(runA)


#### Specify the compartment = B ####
### load intial values
initial_values <-  init.vals()

ibexB <- subset(df, compartment == "B") # subset the data by compartment
str(ibexB)

head(climB)

## Run the simulation ##
w = 1
runB <- run

while(w<10){
runB1 <- gloworm_full(init = initial_values, # initial values
                     ss = ss, ee = ee, lat = latitude, # start and end dates, and location
                     temp = climB$temp2, precip = climB$precip2, # temp and precip parameters
                     species = "teladorsagia", # parameter for parasite species
                     weight = ibexB$totalweight, ages = ibexage$x, 
                     ibexdates = ibexB$date, kgDM = biomass, # parameters for dmi
                     density = ibexB$density, fecs = ibexB$FEC)
runB1$rep <- w
runB <- rbind(runB, runB1)

initial_values$L3faeces <- tail(runB$L3f, 1)
initial_values$L3herbage <- tail(runB$Herbage, 1)
w <- w+1
}

head(runB)



### Specify the compartment = C ####
initial_values <-  init.vals()

ibexC <- subset(df, compartment == "C") # subset the data by compartment
str(ibexC)

head(climC)

## Run the simulation ##
w = 1
runC <- run

while(w<10){
runC1 <- gloworm_full(init = initial_values, # initial values
                     ss = ss, ee = ee, lat = latitude, # start and end dates, and location
                     temp = climC$temp2, precip = climC$precip2, # temp and precip parameters
                     species = "teladorsagia", # parameter for parasite species
                     weight = ibexC$totalweight, ages = ibexage$x, 
                     ibexdates = ibexC$date, kgDM = biomass, # parameters for dmi
                     density = ibexC$density, fecs = ibexC$FEC)
runC1$rep <- w
runC <- rbind(runC, runC1)

initial_values$L3faeces <- tail(runC$L3f, 1)
initial_values$L3herbage <- tail(runC$Herbage, 1)
w <- w+1
}

head(runC)

##### Running the model to plot the effect of ibex ###
### A ibex only
initial_values <-  init.vals()

w = 1
runAibex <- run

while(w<10){
runAibex1 <- gloworm_full(init = initial_values, # initial values
                     ss = ss, ee = ee, lat = latitude, # start and end dates, and location
                     temp = climA$temp2, precip = climA$precip2, # temp and precip parameters
                     species = "teladorsagia", # parameter for parasite species
                     weight = ibexA$weight_ibex, ages = ibexage$x, 
                     ibexdates = ibexA$date, kgDM = biomass, # parameters for dmi
                     density = ibexA$ibexdensity, fecs = ibexA$FEC_ibexnew)
runAibex1$rep <- w
runAibex <- rbind(runAibex, runAibex1)

initial_values$L3faeces <- tail(runAibex$L3f, 1)
initial_values$L3herbage <- tail(runAibex$Herbage, 1)
w <- w+1
}

head(runAibex)


### B ibex only
initial_values <-  init.vals()

w = 1
runBibex <- run

while(w<10){
runBibex1 <- gloworm_full(init = initial_values, # initial values
                     ss = ss, ee = ee, lat = latitude, # start and end dates, and location
                     temp = climB$temp2, precip = climB$precip2, # temp and precip parameters
                     species = "teladorsagia", # parameter for parasite species
                     weight = ibexB$weight_ibex, ages = ibexage$x, 
                     ibexdates = ibexB$date, kgDM = biomass, # parameters for dmi
                     density = ibexB$ibexdensity, fecs = ibexB$FEC_ibexnew)
runBibex1$rep <- w
runBibex <- rbind(runBibex, runBibex1)

initial_values$L3faeces <- tail(runBibex$L3f, 1)
initial_values$L3herbage <- tail(runBibex$Herbage, 1)
w <- w+1
}

head(runBibex)

### C ibex only
initial_values <-  init.vals()

w = 1
runCibex <- run

while(w<10){
runCibex1 <- gloworm_full(init = initial_values, # initial values
                     ss = ss, ee = ee, lat = latitude, # start and end dates, and location
                     temp = climC$temp2, precip = climC$precip2, # temp and precip parameters
                     species = "teladorsagia", # parameter for parasite species
                     weight = ibexC$weight_ibex, ages = ibexage$x, 
                     ibexdates = ibexC$date, kgDM = biomass, # parameters for dmi
                     density = ibexC$ibexdensity, fecs = ibexC$FEC_ibexnew)
runCibex1$rep <- w
runCibex <- rbind(runCibex, runCibex1)

initial_values$L3faeces <- tail(runCibex$L3f, 1)
initial_values$L3herbage <- tail(runCibex$Herbage, 1)
w <- w+1
}

head(runCibex)

###### compile the outputs
runA$compartment <- "A"
runB$compartment <- "B"
runC$compartment <- "C"
runAibex$compartment <- "A"
runBibex$compartment <- "B"
runCibex$compartment <- "C"
runA$Species <- "Alpine ibex and sheep"
runB$Species <- "Alpine ibex and sheep"
runC$Species <- "Alpine ibex and sheep"
runAibex$Species <- "Alpine ibex"
runBibex$Species  <- "Alpine ibex"
runCibex$Species  <- "Alpine ibex"

fullrun <- rbind(runA, runB, runC, runAibex, runBibex, runCibex)

write.csv(fullrun, file = "Modeloutput23_multiple.csv")
fullrun <- read.csv("Modeloutput23_multiple.csv")

head(fullrun)
str(fullrun)

fullrun <- subset(fullrun, rep !=1)

##### Looking at the output ####
fullrun2 <- fullrun %>% 
  dplyr:::group_by(Species) %>% 
  dplyr:::summarize(E = mean(E), L3p = mean(L3p),
                    H = mean(Herbage), Hkg = mean(Herbperkg),  HkgSD = sd(Herbperkg))

##### Looking at the output ####
summary <- fullrun %>% 
  dplyr:::group_by(compartment, Species) %>% 
  dplyr:::summarize(sumE = sum(E), meanE = mean(E), sumL3p = sum(L3p), meanL3p = mean(L3p),
            sumH = sum(Herbage), meanH = mean(Herbage), sdH = sd(Herbage), H = mean(Herbperkg),
            median = median(L3p), 
            quant1 = quantile(L3p, probs = c(.25)),
            quant3 = quantile(L3p, probs = c(.75)))
library(plyr)
library(DescTools)
## AUC
auccalc <- ddply(fullrun, c("compartment", "Species"), function(z){
  AUC(x= z$time, y= z$L3p, method="trapezoid")
})

summary$auc <- auccalc$V1
head(summary)
write.csv(summary, file = "modeloutput_23medianandauc.csv")



#### plot ####
head(fullrun2)
str(fullrun2)
fullrun2$Species <- factor(fullrun2$Species, levels = c("Alpine ibex and sheep", "Alpine ibex"))

e <- ggplot(fullrun2, aes(x = time, y = E, colour = Species))+
  geom_line(size = 0.9, aes(linetype=Species))+
  facet_grid(~compartment)+ ylab("Eggs onto pasture (per ha)")+
  #scale_fill_manual(values=c("lightgrey", "lightblue")) +
  scale_color_manual(values=c("black", "blue2")) +
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+ 
  geom_segment(x = 175, y = -3000, xend = 259, yend = -3000, size = 2, colour = "darkgrey")+
  #geom_segment(aes(x = 259, y = -300, xend = 259, yend = -100),
   #            arrow = arrow(length = unit(0.25, "cm")))+
  theme_classic()+theme(axis.title.x=element_blank(),axis.text.x = element_blank(),
                        legend.position = c(0.15, 0.8))
e

runA$rep <- as.character(runA$rep)

l3 <- ggplot(fullrun2, aes(x = time, y = L3p, colour = Species))+
  geom_line(size = 0.9, aes(linetype=Species))+
  facet_grid(~compartment)+ ylab("L3 on pasture (per ha)")+
  #scale_fill_manual(values=c("lightgrey", "lightblue")) +
  scale_color_manual(values=c("black", "blue2")) +
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+ 
  geom_segment(x = 175, y = -3000, xend = 259, yend = -3000, size = 2, colour = "darkgrey")+
  theme_classic()+theme(axis.title.x=element_blank(),axis.text.x = element_blank(),
                        legend.position = "none")
l3

h <- ggplot(fullrun2, aes(x = time, y = H, colour = Species))+
  geom_line(size = 0.8, aes())+
  facet_grid(~compartment)+ xlab("Month")+ ylab("L3 on herbage (per ha)")+
  #scale_fill_manual(values=c("lightgrey", "lightblue")) +
  scale_color_manual(values=c("black", "blue2")) +
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+ 
  geom_segment(x = 175, y = -300, xend = 259, yend = -300, size = 2, colour = "darkgrey")+
  theme_classic() +theme(axis.text.x = element_text(angle=90, vjust = 0.1),
                         legend.position = "none")
h

png(file="Model output23_multiple.png",width=800,height=1000,res=120)
grid.arrange(e, l3, h, ncol = 1)
dev.off()


#### Grner validation
valid <- read.csv("Grunersimulation_validation.csv")
head(valid)
str(valid)
head(fullrun2)
str(fullrun2)

valid$compartment <- valid$ï..compartment

valid$key <- paste(valid$ï..compartment, valid$time)
fullrun2$key <- paste(fullrun2$compartment, fullrun2$time)
valid$modelherbageboth <- paste(fullrun2[fullrun2$Species == "Alpine ibex and sheep",]$Hkg
                                [match(valid$key, fullrun2[fullrun2$Species == "Alpine ibex and sheep",]$key)])
valid$modelherbageboth <- as.numeric(valid$modelherbageboth)
valid$modelherbagesheep <- valid$modelherbageboth - valid$modelherbageibex


mean(valid$modelherbageboth)
sd(valid$modelherbageboth)

mean(valid$modelherbageibex)
sd(valid$modelherbageibex)

mean(valid$modelherbagesheep)
sd(valid$modelherbagesheep)

mean(valid$herbage)
sd(valid$herbage)

write.csv(valid, file = "Grunersimulation_validation.csv")

herbkgliv_A <- ggplot(fullrun2, aes(x = time, y = Hkg))+
  geom_line(size = 0.6, aes(colour = Species, linetype = Species))+ 
  ylim(-50,450)+
  geom_point(data= valid, aes(y=herbage), size=3)+
  geom_errorbar(data= valid, aes(y=herbage,ymin=herbage-herbsd, ymax=herbage+herbsd))+
  facet_grid(~compartment)+
  scale_color_manual(values=c("black", "darkblue")) +
  scale_fill_manual(values=c("lightgrey", "lightblue")) +
  ylab("L3 per kg herbage")+xlab("Month")+
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                "Sep", "Oct", "Nov", "Dec", "Jan"))+ 
theme_classic()+theme(axis.text.x = element_text(angle=90, vjust = 0.1),
                      legend.position = c(0.85, 0.8))

png(file= "Belledonne_validation_new.png", width = 800, height = 400, res = 120)
herbkgliv_A
dev.off()


##### modeling the validation 
head(valid)

ggplot(valid, aes(y = herbage, x = modelherbageibex))+
  geom_point()+geom_smooth()


mod <- lm(herbage ~ modelherbageboth, data = valid)
summary(mod)

mod2 <- lm(herbage ~ 0 + modelherbageboth, data = valid)
summary(mod2)
confint(mod2, level = 0.95)

mod3 <- lm(herbage ~ modelherbageibex, data = valid)
summary(mod3)

mod4 <- lm(herbage ~ 0 + modelherbageibex, data = valid)
summary(mod4)
confint(mod4, level = 0.95)

cor.test(valid$herbage, valid$modelherbageboth)
cor.test(valid$herbage, valid$modelherbageibex)
cor.test(valid$herbage, valid$modelherbagesheep)


### plot

ggplot(valid, aes(x=modelherbageboth, y=herbage)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_abline(intercept=0, slope=mod2$coefficients[1], color='#2C3E50', size=1.1) + 
  labs(title='Regression through the Origin')

cor.test(valid$modelherbageboth, valid$herbage, method = "pearson")
cor.test(valid$modelherbageibex, valid$herbage, method = "pearson")

valid$modelherbagesheep <- valid$modelherbageboth - valid$modelherbageibex
head(valid)
cor.test(valid$modelherbagesheep, valid$herbage, method = "pearson")


mod5 <- lm(herbage ~ 0 + modelherbagesheep, data = valid)
summary(mod5)
confint(mod5, level = 0.95)

library("ggpubr")
ggscatter(valid, x = "herbage", y = "modelherbageboth", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")
ggscatter(valid, x = "herbage", y = "modelherbageibex", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")
ggscatter(valid, x = "herbage", y = "modelherbagesheep", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")
