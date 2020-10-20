#Last updated 10-9-2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")

#Load required packages------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(lme4)#Linear generated mixed models analysis
library(ochRe)# color palette

#Load Data--------------------------------------------------------------------------
Metadata<- read.csv("Metadata/Metadata-PIPO-NC.csv", na.strings = "N/A", header = TRUE) 

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
            max =max(TC, na.rm = TRUE))%>%
 as.data.frame();TCtsfDes


TNtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TN, na.rm = TRUE), 
            sd = sd(TN, na.rm=TRUE),
            min=min(TN, na.rm = TRUE),
            max =max(TN, na.rm = TRUE))%>%
  as.data.frame();TNtsfDes

TPtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TP, na.rm = TRUE), 
            sd = sd(TP, na.rm=TRUE),
            min=min(TP, na.rm = TRUE),
            max =max(TP, na.rm = TRUE))%>%
  as.data.frame();TPtsfDes


TCTPtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TP, na.rm = TRUE), 
            sd = sd(TC.TP, na.rm=TRUE),
            min=min(TC.TP, na.rm = TRUE),
            max =max(TC.TP, na.rm = TRUE))%>%
  as.data.frame();TCTPtsfDes

TCTNtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TN, na.rm = TRUE), 
            sd = sd(TC.TN, na.rm=TRUE),
            min=min(TC.TN, na.rm = TRUE),
            max =max(TC.TN, na.rm = TRUE))%>%
  as.data.frame();TCTNtsfDes


OMtsfDes<-MetaBurn %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(AvgDepthOM, na.rm = TRUE), 
            sd = sd(AvgDepthOM, na.rm=TRUE),
            min=min(AvgDepthOM, na.rm = TRUE),
            max =max(AvgDepthOM, na.rm = TRUE))%>%
  as.data.frame();OMtsfDes



# EXPORT MEANSTDV TABLE------------------------------------------------------------------------------------
dir.create(file.path("Analysis/Nutrients/Tables/DesStats"), recursive = TRUE)

#Creat one table for all results, but have to add names outside in excel
NutrientsTSF<-cbind(TCtsfDes,TNtsfDes,TPtsfDes,TCTNtsfDes,
                    TCTPtsfDes,OMtsfDes);NutrientsTSF

write.csv(NutrientsTSF, file="Analysis/Nutrients/Tables/DesStats/Nutrient-Des-Stats-TSF.csv") 




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

#Creat one table for all results, but have to add names outside in excel
NutriTSFPerCh<-cbind(TCtsfDes,TNtsfDes,TPtsfDes,TCTNtsfDes,
                    TCTPtsfDes,OMtsfDes);NutriTSFPerCh

write.csv(NutriTSFPerCh, file="Analysis/Nutrients/Tables/PercentChange/Nutrient-PerChange-TSF.csv") 



#
#
#**************************************************************************************----
# --------------------- LINEAR REGRESSION ANALYSIS ----------------------------------------
#**************************************************************************************----
#Test diff poly models to see order of poly to use...if adding model results in 
#non-significance= use previous model

attach(MetaBurn)
X=TimeSinceFire
library(rsq)
library(emmeans)#post hoc

#CARBON---Poly(X, 2) sign ......................................................
fitC<-glm(TC~poly(TimeSinceFire, 2,  raw=TRUE) + (1/Site), data = MetaBurn);summary(fitC)
fitCsum<-summary.glm(fitC);fitCsum
rsq(fitC,adj=TRUE)




#NITROGEN---Poly(X, 2) sign ....................................................
fitN<-glm(TN~poly(X, 2, raw=TRUE)+ (1/Site), data = MetaBurn)
fitNsum<-summary.glm(fitN);fitNsum #signif
rsq(fitN,adj=TRUE)

#PHOSPHORUS---Poly(X, 3) sign ..................................................
fitP<-glm(TP~poly(X, 3, raw=TRUE)+ (1/Site),data=MetaBurn)
fitPsum<-summary.glm(fitP);fitPsum #signif
rsq(fitP,adj=TRUE)

#CARBON TO NITRO---Poly(X, 2) sign .............................................
fitTCTN<-glm(TC.TN ~ poly(X, 2, raw=TRUE)+ (1/Site),data=MetaBurn)

fitTCTNsum<-summary.glm(fitTCTN);fitTCTNsum #signif
rsq(fitTCTN,adj=TRUE)

#CARBON:PHOSPHORUS---Poly(X, 3) sign ...........................................
fitTCTP<-glm(TC.TP~poly(X, 3, raw=TRUE)+ (1/Site), data=MetaBurn)
fitTCTPsum<-summary.glm(fitTCTP);fitTCTPsum#signif
rsq(fitTCTP,adj=TRUE)

#DEPTH OM --- y ~ X sign other lead to nonsignificant intercept.................
fitOM<-glm(AvgDepthOM~X + (1/Site),data=MetaBurn)
fitOMsum<-summary.glm(fitOM);fitOMsum#signif
rsq(fitOM,adj = TRUE)

#eXPORT Restults..............................................................................
dir.create("Analysis/Nutrients/Glm")


capture.output(fitPsum, file = "Analysis/Nutrients/Glm/Carbon-Regression.csv")
capture.output(fitNsum, file = "Analysis/Nutrients/Glm/Nitrogen-Regression.csv")
capture.output(fitPsum, file = "Analysis/Nutrients/Glm/Phosphorus-Regression.csv")
capture.output(fitTCTNsum, file = "Analysis/Nutrients/Glm/TC.TN-Regression.csv")
capture.output(fitTCTPsum, file = "Analysis/Nutrients/Glm/TC.TP-Regression.csv")
capture.output(fitOMsum, file = "Analysis/Nutrients/Glm/AvgDepthOM-Regression.csv")




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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,angle=90,hjust=0.5),
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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,angle=90,hjust=0.5),
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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=18,angle=90,hjust=0.5),
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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,angle=90,hjust=0.5),
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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,angle=90,hjust=0.5),
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
        axis.title.y = element_text(margin = margin(t=0,r=4,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=0,l=0)),
        text = element_text(size=23),
        axis.text.y = element_text(size=19, angle=90,hjust=0.5),
        axis.text.x = element_text(size=18),
        legend.position="none")+ylim(0,6)+
  labs(x = "", y= "Depth of OM (cm)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11));glmOM
  


dir.create(file.path("Analysis/Nutrients/Graphs/"), recursive = TRUE)
pdf("Analysis/Nutrients/Graphs/LinearRegressionsTimeSinceFire.pdf", height=8.5, width=14,onefile=FALSE)
ggarrange(glmTP,glmOM,glmTCTN,glmTCTP,glmTC,glmTN,ncol=3,nrow=2, common.legend = TRUE, 
          legend="bottom", labels = "auto", align="hv", label.x=0, label.y=1,
          font.label = list(size=18,face="bold"))
dev.off()













