
#RESET R.............................................................
rm(list=ls())

#Set working directory------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

#Load required packages----------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(lme4)#Linear generated mixed models analysis
library(ochRe)# color palette

#Load Data--------------------------------------------------------------------------
Metadata<- read.csv("Metadata-PIPO-NC.csv", na.strings = "N/A", header = TRUE) 



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
#

attach(Metadata)

#TSF..............................................................
TCDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC, na.rm = TRUE), 
            sd = sd(TC, na.rm=TRUE),
            min=min(TC, na.rm = TRUE),
            max =max(TC, na.rm = TRUE));TCDesTrt


TNDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TN, na.rm = TRUE), 
            sd = sd(TN, na.rm=TRUE),
            min=min(TN, na.rm = TRUE),
            max =max(TN, na.rm = TRUE));TNDesTrt

TPDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TP, na.rm = TRUE), 
            sd = sd(TP, na.rm=TRUE),
            min=min(TP, na.rm = TRUE),
            max =max(TP, na.rm = TRUE));TPDesTrt


TCTPDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TP, na.rm = TRUE), 
            sd = sd(TC.TP, na.rm=TRUE),
            min=min(TC.TP, na.rm = TRUE),
            max =max(TC.TP, na.rm = TRUE));TCTPDesTrt

TCTNDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(TC.TN, na.rm = TRUE), 
            sd = sd(TC.TN, na.rm=TRUE),
            min=min(TC.TN, na.rm = TRUE),
            max =max(TC.TN, na.rm = TRUE));TCTNDesTrt


OMDesTrt<-Metadata %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(AvgDepthOM, na.rm = TRUE), 
            sd = sd(AvgDepthOM, na.rm=TRUE),
            min=min(AvgDepthOM, na.rm = TRUE),
            max =max(AvgDepthOM, na.rm = TRUE));OMDesTrt



# EXPORT MEANSTDV TABLE------------------------------------------------------------------------------------
dir.create("Analysis/Nutrients/Tables/DesStats")

write.csv(TCDesTrt, file="Analysis/Nutrients/Tables/DesStats/Carbon-Des-Stats-TRT.csv") 
write.csv(TNDesTrt, file="Analysis/Nutrients/Tables/DesStats/Nitrogen-Des-Stats-TRT.csv") 
write.csv(TPDesTrt, file="Analysis/Nutrients/Tables/DesStats/Phosphorus-Des-Stats-TRT.csv") 
write.csv(TCTNDesTrt, file="Analysis/Nutrients/Tables/DesStats/CarbonNitrogen-Des-Stats-TRT.csv") 
write.csv(TCTPDesTrt, file="Analysis/Nutrients/Tables/DesStats/CarbonPhosphorus-Des-Stats-TRT.csv") 
write.csv(OMDesTrt, file="Analysis/Nutrients/Tables/DesStats/OM-Des-Stats-TRT.csv") 


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

write.csv(TCPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/Carbon-Per-Change-TRT.csv") 
write.csv(TNPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/Nitrogen-Per-Change-TRT.csv") 
write.csv(TPPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/Phosphorus-Per-Change-TRT.csv") 
write.csv(TCTNPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/CarbonNitrogen-Per-Change-TRT.csv") 
write.csv(TCTPPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/CarbonPhosphorus-Per-Change-TRT.csv") 
write.csv(OMPerChgTrt, file="Analysis/Nutrients/Tables/PercentChange/OM-Per-Change-TRT.csv") 

#
#
#***************************************************************************************************************----
# ---------------- TEST FOR SIGNIFICANCE ---------------------------------------------------------------------------
#***************************************************************************************************************----

#KRUSKAL WALLIS TEST AND POST-HOC WILCOXON ...................................................

#CARBON.................................................................
KWcarbon<-kruskal.test(TC ~ Treatment);KWcarbon #p=0.4344

KWnitrogen<-kruskal.test(TN ~ Treatment);KWnitrogen#p=0.1696

KWphosphorus<-kruskal.test(TP ~ Treatment);KWphosphorus #p=3.98e-14

KWtctn<-kruskal.test(TC.TN ~ Treatment);KWtctn #p=0.02351

KWtctp<-kruskal.test(TC.TP ~ Treatment);KWtctp #p=1.356e-06

KWom<-kruskal.test(AvgDepthOM ~ Treatment);KWom #p<2.2e-16



#
#
#***************************************************************************************************************----
# ---------------- CREATE PLOTS ------------------------------------------------------------------------------------
#***************************************************************************************************************----
#

TcTrt<-ggplot(Metadata, aes(x=Treatment, y=TC)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,10)+  
  labs(x = "", y= "Total Carbon (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TcTrt
TcTrt1<-TcTrt + stat_compare_means(method = "kruskal", label.y = 10,); TcTrt1



TnTrt<-ggplot(Metadata, aes(x=Treatment, y=TN)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,2.5)+  
  labs(x = "", y= "Total Nitrogen (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TnTrt
TnTrt1<-TnTrt + stat_compare_means(method = "kruskal", label.y = 2.5); TnTrt1


TpTrt<-ggplot(Metadata, aes(x=Treatment, y=TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,0.3)+  
  labs(x = "", y= "Total Phosphorus (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TpTrt
TpTrt1<-TpTrt + stat_compare_means(method = "kruskal", label.y = 0.3); TpTrt1



TcTpTrt<-ggplot(Metadata, aes(x=Treatment, y=TC.TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=16,hjust=0.5),
                   axis.text.x = element_text(size=16))+ylim(0,250)+
  labs(x = "", y= "TC:TP Ratio (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); TcTpTrt
TcTpTrt1<-TcTpTrt + stat_compare_means(method = "kruskal", label.y = 250);TcTpTrt1


TcTnTrt<-ggplot(Metadata, aes(x=Treatment, y=TC.TP)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+  ylim(0,250)+
  labs(x = "", y= "TC:TN Ratio (%)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));TcTnTrt
TcTnTrt1<-TcTnTrt + stat_compare_means(method = "kruskal", label.y = 250); TcTnTrt1



OMTrt<-ggplot(Metadata, aes(x=Treatment, y=AvgDepthOM)) +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ 
  labs(x = "", y= "Avg Depth OM (cm)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));OMTrt
OMTrt1<-OMTrt + stat_compare_means(method = "kruskal", label.y = 15);OMTrt1




#EXPORT GRAPHS-------------------------------------------------------------------------------------

pdf("Analysis/Nutrients/Graphs/Nutrients-Treatment.pdf", height=9, width=12,onefile=FALSE)
ggarrange(TpTrt1,OMTrt1,TcTpTrt1,TcTnTrt1,TcTrt1,TnTrt1,ncol=3,nrow=2, common.legend = TRUE, 
          legend="bottom", labels = "auto", align="hv",
          font.label = list(size=16,face="bold"))
dev.off()








attach(Metadata)
plot(TC~Treatment, ylab="Total Carbon (%)", col=c("slategray", "burlywood"))
plot(TN~Treatment, ylab= "Total Nitrogen (%)",col=c("slategray", "burlywood"))
plot(TP~Treatment, ylab="Total Phosphorus (%)", col=c("slategray", "burlywood"))
plot(TC.TP~Treatment, ylab= "Ratio TC:TN (%)",col=c("slategray", "burlywood"))
plot(TC.TN~Treatment, ylab="Ratio TC:TP (%)", col=c("slategray", "burlywood"))
plot(AvgDepthOM~Treatment, ylab="Average depth OM (cm)", col=c("slategray", "burlywood"))



