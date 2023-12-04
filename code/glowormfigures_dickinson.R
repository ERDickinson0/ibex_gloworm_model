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
