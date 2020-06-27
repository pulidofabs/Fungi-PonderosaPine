#Last updated June 23, 2020


#Set working directory-------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

#Load required packages------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(lme4)#Linear generated mixed models analysis
library(ochRe)# color palette

#Load Data--------------------------------------------------------------------------
Metadata<- read.csv("Metadata-PIPO-NC.csv", na.strings = "N/A", header = TRUE) 

#Quality control- ------------------------------------------------------------------
#Subset metadata to maintain only burned samples as we are interested in the 
#changes in the nutrients in the burned plots 

attach(Metadata)
MetaBurn<-Metadata[which(Treatment== "Burned"), ]
head(MetaBurn)
dim(MetaBurn)#107x26
str(MetaBurn)
detach(Metadata)


#----
#----
#*************************************************************************************************----
# --------------------   Descriptive statistics ------------------------------------------------------
#*************************************************************************************************----
attach(MetaBurn)

#TSF..............................................................
TCtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TC, na.rm = TRUE), 
            sd = sd(TC, na.rm=TRUE),
            min=min(TC, na.rm = TRUE),
            max =max(TC, na.rm = TRUE));TCtsfDes


TNtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TN, na.rm = TRUE), 
            sd = sd(TN, na.rm=TRUE),
            min=min(TN, na.rm = TRUE),
            max =max(TN, na.rm = TRUE));TNtsfDes

TPtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TP, na.rm = TRUE), 
            sd = sd(TP, na.rm=TRUE),
            min=min(TP, na.rm = TRUE),
            max =max(TP, na.rm = TRUE));TPtsfDes


TCTPtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TP, na.rm = TRUE), 
            sd = sd(TC.TP, na.rm=TRUE),
            min=min(TC.TP, na.rm = TRUE),
            max =max(TC.TP, na.rm = TRUE));TCTPtsfDes

TCTNtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TN, na.rm = TRUE), 
            sd = sd(TC.TN, na.rm=TRUE),
            min=min(TC.TN, na.rm = TRUE),
            max =max(TC.TN, na.rm = TRUE));TCTNtsfDes


OMtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(AvgDepthOM, na.rm = TRUE), 
            sd = sd(AvgDepthOM, na.rm=TRUE),
            min=min(AvgDepthOM, na.rm = TRUE),
            max =max(AvgDepthOM, na.rm = TRUE));OMtsfDes



# EXPORT MEANSTDV TABLE------------------------------------------------------------------------------------
dir.create("Analysis/Nutrients/Tables/DesStats")

write.csv(TCtsfDes, file="Analysis/Nutrients/Tables/DesStats/CarbonDes-Stats-TSF.csv") 
write.csv(TNtsfDes, file="Analysis/Nutrients/Tables/DesStats/NitrogenDes-Stats-TSF.csv") 
write.csv(TPtsfDes, file="Analysis/Nutrients/Tables/DesStats/PhosphorusDes-Stats-TSF.csv") 
write.csv(TCTNtsfDes, file="Analysis/Nutrients/Tables/DesStats/CarbonNitrogenDes-Stats-TSF.csv") 
write.csv(TCTPtsfDes, file="Analysis/Nutrients/Tables/DesStats/CarbonPhosphorusDes-Stats-TSF.csv") 
write.csv(OMtsfDes, file="Analysis/Nutrients/Tables/DesStats/OMDes-Stats-TSF.csv") 


#
#
#***********************************************************************************************************************----
#--- PERCENT CHANGE CALCULATIONS -------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
library(dplyr)

TCPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(TC,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCPerChg

TNPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(TN,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TNPerChg

TPPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(TP,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TPPerChg

TCTPPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(TC.TP,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCTPPerChg

TCTNPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(TC.TN,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCTNPerChg

OMPerChg<-MetaBurn%>%
  group_by(TimeSinceFire) %>%
  summarise(mean = mean(AvgDepthOM,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();OMPerChg


#EXPORT TABLES -----------------------------------------------------------------------
dir.create("Analysis/Nutrients/Tables/PercentChange")

write.csv(TCPerChg, file="Analysis/Nutrients/Tables/PercentChange/Carbon-Per-Change-TSF.csv") 
write.csv(TNPerChg, file="Analysis/Nutrients/Tables/PercentChange/Nitrogen-Per-Change-TSF.csv") 
write.csv(TPPerChg, file="Analysis/Nutrients/Tables/PercentChange/Phosphorus-Per-Change-TSF.csv") 
write.csv(TCTNPerChg, file="Analysis/Nutrients/Tables/PercentChange/CarbonNitrogen-Per-Change-TSF.csv") 
write.csv(TCTPPerChg, file="Analysis/Nutrients/Tables/PercentChange/CarbonPhosphorus-Per-Change-TSF.csv") 
write.csv(OMPerChg, file="Analysis/Nutrients/Tables/PercentChange/OM-Per-Change-TSF.csv") 



#
#
#**************************************************************************************----
# --------------------- LINEAR REGRESSION ANALYSIS ----------------------------------------
#**************************************************************************************----
#Test diff poly models to see order of poly to use...if adding model results in 
#non-significance= use previous model

attach(MetaBurn)
X=TimeSinceFire


#CARBON---Poly(X, 2) sign ......................................................
fitC<-glm(TC~poly(X, 2, raw=TRUE), data = MetaBurn )
summary(fitC)


#NITROGEN---Poly(X, 2) sign ....................................................
fitN<-glm(TN~poly(X, 2, raw=TRUE), data = MetaBurn)
summary(fitN)#*


#PHOSPHORUS---Poly(X, 3) sign ..................................................
fitP<-glm(TP~poly(X, 3, raw=TRUE),data=MetaBurn)
summary(fitP) #sign


#CARBON:PHOSPHORUS---Poly(X, 3) sign ...........................................
fitTCTP<-glm(TC.TP~poly(X, 3, raw=TRUE), data=MetaBurn)
summary(fitTCTP)#sig


#CARBON TO NITRO---Poly(X, 2) sign .............................................
fitTCTN<-glm(TC.TN ~ poly(X, 2, raw=TRUE),data=MetaBurn)
summary(fitTCTN)#sig


#DEPTH OM --- y ~ X sign other lead to nonsignificant intercept.................
fitOM<-glm(AvgDepthOM~X, data=MetaBurn)
summary(fitOM)#sig



#eXPORT Restults..............................................................................
dir.create("Analysis/Nutrients/Glm")

capture.output(summary(fitP), file = "Analysis/Nutrients/Glm/Carbon-Regression.csv")
capture.output(summary(fitN), file = "Analysis/Nutrients/Glm/Nitrogen-Regression.csv")
capture.output(summary(fitP), file = "Analysis/Nutrients/Glm/Phosphorus-Regression.csv")
capture.output(summary(fitTCTN), file = "Analysis/Nutrients/Glm/TC.TN-Regression.csv")
capture.output(summary(fitTCTP), file = "Analysis/Nutrients/Glm/TC.TP-Regression.csv")
capture.output(summary(fitOM), file = "Analysis/Nutrients/Glm/AvgDepthOM-Regression.csv")




#
#
#***************************************************************************************************----
#------------ PLOT RESULTS------------------------------------------------------------------------------
#***************************************************************************************************----
(mfrow=c(1,1))
library(ochRe)


#CARBON......................................................................................----
formula1 <- y ~ poly(x, 2, raw = TRUE)

glmTC<-ggplot(MetaBurn, aes(TimeSinceFire, TC))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula1, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
        formula = formula1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="")+ ylim(0,10)+
  labs(x = "", y= "Total Carbon (%)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmTC


#NITROGEN...................................................................................----
formula2 <- y~poly(x, 2, raw = TRUE)

glmTN<-ggplot(MetaBurn, aes(TimeSinceFire, TN))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula2, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
        formula = formula2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ ylim(-0.1,3)+
  labs(x = "Time Since Fire (years)", y= "Total Nitrogen (%)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmTN
  


#PHOSPHORUS...............................................................................----
formula3<- y ~ poly(x,3, raw = TRUE)

glmTP<-ggplot(MetaBurn, aes(TimeSinceFire, TP))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula3, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
         formula = formula3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=18,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ ylim(0,0.4)+
  labs(x = "", y= "Total Phosphorus (%)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmTP

  

#CARBON+PHOSPHORUS.........................................................................----
formula4<-y~poly(x,3, raw = TRUE)

glmTCTP<-ggplot(MetaBurn, aes(TimeSinceFire, TC.TP))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula4, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
        formula = formula4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ ylim(-90,130)+
  labs(x = "Time Since Fire (years)", y= "TC:TP Ratio (%)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmTCTP

  
#CARBON+NITROGEN.........................................................................----
formula5 <- y ~ poly(x, 2, raw = TRUE)

glmTCTN<-ggplot(MetaBurn, aes(TimeSinceFire, TC.TN))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula5, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
        formula = formula5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ ylim(0,30)+
  labs(x = "", y= "TC:TN Ratio (%)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmTCTN
  
  
#ORGANIC MATTER DEPTH .........................................................................----

formula6<- y~poly(x, raw = TRUE)

glmOM<-ggplot(MetaBurn, aes(TimeSinceFire, AvgDepthOM))+
  geom_point(aes(col= Site), size=1.5)+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method='glm', formula=formula6, se=TRUE, col="blue", lwd=1.5)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
        formula = formula6)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ylim(0,6)+
  labs(x = "", y= "Depth of OM (cm)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmOM
  

pdf("Analysis/Nutrients/Graphs/LinearRegressionsTimeSinceFire.pdf", height=8.5, width=12,onefile=FALSE)
ggarrange(glmTP,glmOM,glmTCTP,glmTCTN,glmTC,glmTN, ncol=3,nrow=2, common.legend = TRUE, 
          legend="bottom", labels = "auto", align="hv", label.x=0, label.y=1,
          font.label = list(size=18,face="bold"))
dev.off()













