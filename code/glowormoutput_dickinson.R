###############################################################################
#                         Gloworm simulations                                 #
#                                                                             #
#                       Author: Eleanor Dickinson                             #
#                       Date: November 2023                                   #
#                                                                             #
# 1. Summarizing model output.                                                #
#                                                                             #
# 3. Calculating host exposure.                                               #
#                                                                             #
# 2. Validation code comparing the model output to pasture counts reported by #
#  Gruner et al., 2008.                                                       #
###############################################################################

library(plyr)
library(dplyr)
library(DescTools)

fullrun <- read.csv("Modeloutput.csv")

###############################################################################
# Summarizing model output                                                    #
###############################################################################

##### Summarizing ####
summary <- fullrun %>% 
  dplyr:::group_by(compartment, Species) %>% 
  dplyr:::summarize(meanE = mean(E), meanL3p = mean(L3p),
                    meanH = mean(Herbage), sdH = sd(Herbage), 
                    median = median(L3p, na.rm = T),
                    QR1 = quantile(L3p, c(0.25), na.rm = T),
                    QR3 = quantile(L3p, c(0.75), na.rm = T))

##### Calculating AUC ####
auccalc <- ddply(fullrun, c("compartment", "Species"), function(z){
  AUC(x= z$time, y= z$L3p, method="trapezoid")
})
summary$auc <- auccalc$V1
summary
#------------------------------------------------------------------------------

###############################################################################
# Calculating host exposure                                                   #
###############################################################################

#### Total AUC in the period and compartment occupied by each host. 
## Alpine ibex 
ibex <- subset(fullrun, Species == "Alpine ibex and sheep")

ibexa <- subset(ibex, compartment == "A")
ibexa <- subset(ibexa, time <= 161 | time >= 245)

ibexb <- subset(ibex, compartment == "B")

ibexc <- subset(ibex, compartment == "C")
ibexc1 <- subset(ibexc, time >= 161 & time <= 308)
ibexc2 <- subset(ibexc, time >= 322 & time <= 329)
ibexc3 <- subset(ibexc, time >= 336)
ibexc <- rbind(ibexc1, ibexc2, ibexc3)

ibexf <- rbind(ibexa, ibexb, ibexc)

summaryibex <- ibexf %>% 
  dplyr:::group_by(compartment) %>% 
  dplyr:::summarize(median = median(L3p, na.rm = T),
                    QR1 = quantile(L3p, c(0.25), na.rm = T),
                    QR3 = quantile(L3p, c(0.75), na.rm = T))
auccalc <- ddply(ibexf, c("compartment", "Species"), function(z){
  AUC(x= z$time, y= z$L3p, method="trapezoid")
})
summaryibex$auc <- auccalc$V1
summaryibex$Species <- "Ibex"

## Sheep 
ibex <- subset(fullrun, Species == "Alpine ibex and sheep")

ibexa <- subset(ibex, compartment == "A")
ibexa <- subset(ibexa, time >= 182 & time <= 252)

ibexb <- subset(ibex, compartment == "B")
ibexb <- subset(ibexb, time >= 182 & time <= 252)

ibexc <- subset(ibex, compartment == "C")
ibexc <- subset(ibexc, time >= 196 & time <= 252)

ibexf <- rbind(ibexa, ibexb, ibexc)

summarysheep <- ibexf %>% 
  dplyr:::group_by(compartment) %>% 
  dplyr:::summarize(median = median(L3p, na.rm = T),
                    QR1 = quantile(L3p, c(0.25), na.rm = T),
                    QR3 = quantile(L3p, c(0.75), na.rm = T))
auccalc <- ddply(ibexf, c("compartment", "Species"), function(z){
  AUC(x= z$time, y= z$L3p, method="trapezoid")
})
summarysheep$auc <- auccalc$V1
summarysheep$Species <- "Sheep"

summaryfull <- rbind(summaryibex, summarysheep)
summaryfull
#------------------------------------------------------------------------------

###############################################################################
#  Validation code comparing the model output to pasture counts reported by   #
#  Gruner et al., 2008.                                                       #
###############################################################################

valid <- read.csv("Gruner_validation.csv")
fullrun <- read.csv("Modeloutput.csv")

##### Summarising the output ####
fullrun2 <- fullrun %>% 
  dplyr:::group_by(Species) %>% 
  dplyr:::summarize(E = mean(E), L3p = mean(L3p),
                    H = mean(Herbage), Hkg = mean(Herbperkg),  
                    HkgSD = sd(Herbperkg))

### create a key to paste and match the data 
valid$key <- paste(valid$compartment, valid$time)
fullrun2$key <- paste(fullrun2$compartment, fullrun2$time)

### paste and match the data
valid$modelherbageboth <- paste(fullrun2[fullrun2$Species == 
                                           "Alpine ibex and sheep",]$Hkg
                                [match(valid$key, fullrun2[fullrun2$Species == 
                                         "Alpine ibex and sheep",]$key)])
valid$modelherbageboth <- as.numeric(valid$modelherbageboth)

### calculate the number of L3 on herbage from sheep
valid$modelherbagesheep <- valid$modelherbageboth - valid$modelherbageibex


#### Testing the validation ####
cor.test(valid$herbage, valid$modelherbageboth, method = "pearson")
cor.test(valid$herbage, valid$modelherbageibex, method = "pearson")
cor.test(valid$herbage, valid$modelherbagesheep, method = "pearson")


#### Plot the vaidation over the model output ####
herbkgliv_A <- ggplot(fullrun2, aes(x = time, y = Hkg))+
  geom_line(size = 0.6, aes(colour = Species, linetype = Species))+ 
  ylim(-50,450)+
  geom_point(data= valid, aes(y=herbage), size=3)+
  geom_errorbar(data= valid, aes(y=herbage,ymin=herbage-herbsd, 
                                 ymax=herbage+herbsd))+
  facet_grid(~compartment)+
  scale_color_manual(values=c("black", "darkblue")) +
  scale_fill_manual(values=c("lightgrey", "lightblue")) +
  ylab("L3 per kg herbage")+xlab("Month")+
  scale_x_continuous(breaks = c(0,30,60,90,120,150,180,210,240,270,300,330,360),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", 
                                "Jan"))+ 
  theme_classic()+theme(axis.text.x = element_text(angle=90, vjust = 0.1),
                        legend.position = c(0.85, 0.8))

herbkgliv_A

#------------------------------------------------------------------------------
