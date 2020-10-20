#10-9-2020

#RESET R.............................................................
rm(list=ls())

#Set working directory------------------------------------
setwd("C:/Users/juchu/Dropbox/6-PIPO")
setwd("C:/Users/fabipc/Dropbox/6-PIPO")

#Load required packages----------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(lme4)#Linear generated mixed models analysis
library(ochRe)# color palette

#Load Data--------------------------------------------------------------------------
Metadata<- read.csv("Metadata/Metadata-PIPO-NC.csv", na.strings = "N/A", header = TRUE) 
#
#
#***************************************************************************************************************----
# ---------------- QUALITY CONTROL ---------------------------------------------------------------------------
#***************************************************************************************************************----
#Set reference level to unburned....................................................
Metadata$Treatment <- try(relevel(Metadata$Treatment , "Unburned"))
levels(Metadata$Treatment)


#NORMALITY .........................................................................
#Data cannot be transformed only one but decided to proceed with raw data
#Data not normal
attach(Metadata)

# Transformations of data.......................
TC1<-log(TC)
TN1<-log(TN)
CN1<-(TC.TN^3)
TP1<-(TP)^3
OM1<-sqrt(AvgDepthOM)
CP1<-log(TC.TP+1)

#test for normality.............................
shapiro.test(TC1)#normal
shapiro.test(TN1)#NOT 
shapiro.test(CN1)#NOT
shapiro.test(TP1)#NOT
shapiro.test(OM1)#NOT
shapiro.test(CP1)#NOT


#
#
#***************************************************************************************************************----
# ---------------- DESCRIPTIVE STATISTICS --------------------------------------------------------------------------
#***************************************************************************************************************----
attach(Metadata)

#TSF..............................................................
TCDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC, na.rm = TRUE), 
            sd = sd(TC, na.rm=TRUE),
            min=min(TC, na.rm = TRUE),
            max =max(TC, na.rm = TRUE))%>%
            as.data.frame();TCDesTrt


TNDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TN, na.rm = TRUE), 
            sd = sd(TN, na.rm=TRUE),
            min=min(TN, na.rm = TRUE),
            max =max(TN, na.rm = TRUE))%>%
            as.data.frame();TNDesTrt

TPDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TP, na.rm = TRUE), 
            sd = sd(TP, na.rm=TRUE),
            min=min(TP, na.rm = TRUE),
            max =max(TP, na.rm = TRUE))%>%
           as.data.frame();TPDesTrt


TCTPDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TP, na.rm = TRUE), 
            sd = sd(TC.TP, na.rm=TRUE),
            min=min(TC.TP, na.rm = TRUE),
            max =max(TC.TP, na.rm = TRUE))%>%
            as.data.frame();TCTPDesTrt

TCTNDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TN, na.rm = TRUE), 
            sd = sd(TC.TN, na.rm=TRUE),
            min=min(TC.TN, na.rm = TRUE),
            max =max(TC.TN, na.rm = TRUE))%>%
          as.data.frame();TCTNDesTrt


OMDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(AvgDepthOM, na.rm = TRUE), 
            sd = sd(AvgDepthOM, na.rm=TRUE),
            min=min(AvgDepthOM, na.rm = TRUE),
            max =max(AvgDepthOM, na.rm = TRUE))%>%
            as.data.frame();OMDesTrt

#Creat one table for all results, but have to add names outside in excel
EV<-cbind(TCDesTrt,TNDesTrt,TPDesTrt,TCTNDesTrt,TCTPDesTrt,OMDesTrt);EV

# EXPORT MEANSTDV TABLE------------------------------------------------------------------------------------
dir.create("Analysis/Nutrients/Tables/DesStats")

write.csv(EV, file="Analysis/Nutrients/Tables/DesStats/Nutrient-Des-Stats-TRT.csv") 


#
#
#***************************************************************************************************************----
# ---------------- PERCENT CHANGE CALCULATIONS ---------------------------------------------------------------------
#***************************************************************************************************************----
#

library(dplyr)

TCPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(TC,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCPerChgTrt

TNPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(TN,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TNPerChgTrt

TPPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(TP,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TPPerChgTrt

TCTPPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(TC.TP,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCTPPerChgTrt

TCTNPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(TC.TN,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TCTNPerChgTrt

OMPerChgTrt<-Metadata%>%
  group_by(Treatment) %>%
  summarise(mean = mean(AvgDepthOM,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();OMPerChgTrt


#EXPORT TABLES -----------------------------------------------------------------------
dir.create("Analysis/Nutrients/Tables/PercentChange")

#Creat one table for all results, but have to add names outside in excel
PerChange<-cbind(TCPerChgTrt,TNPerChgTrt,TPPerChgTrt,TCTNPerChgTrt,
                 TCTPPerChgTrt,OMPerChgTrt);PerChange

write.csv(PerChange,"Analysis/Nutrients/Tables/PercentChange/Nutrient-PercentChange.csv")

#
#
#***************************************************************************************************************----
# ---------------- TEST FOR SIGNIFICANCE ---------------------------------------------------------------------------
#***************************************************************************************************************----

attach(Metadata)
glm1<- glm(TC ~ Treatment +(1/Site/Plot/Subplot), data = Metadata);summary(glm1)
glm2<- glm(TC ~ Treatment + (1/Site/Plot), data = Metadata);summary(glm2)
glm3<- glm(TC ~ Treatment + (1/Site), data = Metadata);summary(glm3)#no conver
glm4<- glm(TC ~ Treatment + (1/Plot), data = Metadata);summary(glm4)

AIC(glm1,glm2)#same aic
AIC(glm2,glm3)
AIC(glm3,glm4)

anova(glm1,glm2)

#Will do a univariate regression using glm since I can use non normal data
#GLM models dependent variable Yi does NOT need to be normally distributed
glmTC<- glm(TC ~ Treatment + (1/Site), data = Metadata);summary(glmTC)#0.658
glmTN<- glm(TN ~ Treatment + (1/Site), data = Metadata);summary(glmTN)#0.00804 **
glmTP<- glm(TP ~ Treatment + (1/Site), data = Metadata);summary(glmTP)#8.21e-15 ***
glmTCTP<- glm(TC.TP ~ Treatment + (1/Site), data = Metadata);summary(glmTCTP)#0.000118 ***
glmTCTN<- glm(TC.TN ~ Treatment + (1/Site), data = Metadata);summary(glmTCTN)#0.000666 ***
glmOM<- glm(AvgDepthOM ~ Treatment + (1/Site), data = Metadata);summary(glmOM)#<2e-16


library(Publish)
publish(glmTC)


#Export GLM results -------------------------------------------------------------------------
capture.output(summary(glmTC), file="Analysis/Nutrients/Tables/TC-glm-Trt.csv")
capture.output(summary(glmTN), file="Analysis/Nutrients/Tables/TN-glm-Trt.csv")
capture.output(summary(glmTP), file="Analysis/Nutrients/Tables/TP-glm-Trt.csv")
capture.output(summary(glmTCTN), file="Analysis/Nutrients/Tables/TCTN-glm-Trt.csv")
capture.output(summary(glmTCTP), file="Analysis/Nutrients/Tables/TCTP-glm-Trt.csv")
capture.output(summary(glmOM), file="Analysis/Nutrients/Tables/OM-glm-Trt.csv")




#
#
#***************************************************************************************************************----
# ---------------- CREATE PLOTS ------------------------------------------------------------------------------------
#***************************************************************************************************************----
#

TcTrt<-ggplot(Metadata, aes(x=Treatment, y=TC)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ ylim(0,10)+  
  labs(x = "", y= "Total Carbon (%)")+ 
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TcTrt
TcTrt1<-TcTrt + stat_compare_means(method = "anova", label.y = 10,); TcTrt1



TnTrt<-ggplot(Metadata, aes(x=Treatment, y=TN)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,2.5)+  
  labs(x = "", y= "Total Nitrogen (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TnTrt
TnTrt1<-TnTrt + stat_compare_means(method = "anova", label.y = 2.5); TnTrt1


TpTrt<-ggplot(Metadata, aes(x=Treatment, y=TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,0.3)+  
  labs(x = "", y= "Total Phosphorus (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TpTrt
TpTrt1<-TpTrt + stat_compare_means(method = "anova", label.y = 0.3); TpTrt1



TcTpTrt<-ggplot(Metadata, aes(x=Treatment, y=TC.TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
         panel.border = element_rect(colour = "black", fill=NA, size=.7), 
         text = element_text(size=23),
         axis.text.y = element_text(size=16,hjust=0.5),
         axis.text.x = element_text(size=16))+ylim(0,250)+
  labs(x = "", y= "TC:TP Ratio (%)")+ 
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TcTpTrt
TcTpTrt1<-TcTpTrt + stat_compare_means(method = "anova", label.y = 250);TcTpTrt1


TcTnTrt<-ggplot(Metadata, aes(x=Treatment, y=TC.TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+  ylim(0,250)+
  labs(x = "", y= "TC:TN Ratio (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));TcTnTrt
TcTnTrt1<-TcTnTrt + stat_compare_means(method = "anova", label.y = 250); TcTnTrt1



OMTrt<-ggplot(Metadata, aes(x=Treatment, y=AvgDepthOM)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  labs(x = "", y= "Avg Depth OM (cm)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));OMTrt
OMTrt1<-OMTrt + stat_compare_means(method = "anova", label.y = 15);OMTrt1




#EXPORT GRAPHS-------------------------------------------------------------------------------------

pdf("Analysis/Nutrients/Graphs/Nutrients-Treatment-new.pdf", height=9, width=12,onefile=FALSE)
ggarrange(TpTrt1,TnTrt1,TcTnTrt1,TcTpTrt1,OMTrt1,TcTrt1,
          ncol=3,nrow=2, common.legend = TRUE, 
          legend="bottom", labels = "auto", align="hv",
          font.label = list(size=16,face="bold"))
dev.off()




