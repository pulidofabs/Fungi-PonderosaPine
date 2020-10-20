#Sept 15, 2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY--------------------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
library(ggplot2)
library(ggpubr)
library(dplyr)


#Load Data--------------------------------------------------------------------------------------------------------
MetaRareSap<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobeNew.csv", na.strings = "N/A", header = TRUE) 
SaprobeRare3<-read.csv("Analysis/ASV-Tables/Saprobes/RareSaprobe.csv", row.names = 1, check.names = FALSE)


##
##
#*****************************************************************************************************************----
# QUALITY CONTROL )---------------------------------------------------------------------------------------------------
#*****************************************************************************************************************----
# Convert to factor...........................................................
MetaRareSap$FireYear <- as.factor(MetaRareSap$FireYear)
MetaRareSap$TimeSinceFire <- as.factor(MetaRareSap$TimeSinceFire)
MetaRareSap$Treatment<-as.factor(MetaRareSap$Treatment)


# Set reference level of data for comparisons ..............................................................
levels(MetaRareSap$Treatment)
MetaRareSap$Treatment <- try(relevel(MetaRareSap$Treatment , "Unburned"))
levels(MetaRareSap$Treatment) # in correct order


# SUBSET DATA ----------------------------------------------------------------------------------------------------
# * * UnburnedSap PLOTS ............................................
attach(MetaRareSap)
UnburnedSap<-MetaRareSap[which(Treatment== "Unburned"), ]
head(UnburnedSap[1:2,]);dim(UnburnedSap)#102x34

# * * BurnedSap PLOTS ..............................................
BurnedSap<-MetaRareSap[which(Treatment== "Burned"), ]
head(BurnedSap[1:2,]);dim(BurnedSap)#106x34
detach(MetaRareSap)


#
#*****************************************************************************************************************************----
#DESRIPTIVE STATISTICS----------------------------------------------------------------------------------------------------------
#*****************************************************************************************************************************----

#Treatment and timesincefire on all taxa .........................................................
TrtStats<-MetaRareSap %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE),
            sd = sd(S.obs), 
            se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TrtStats # B=8.413793, U=15.514286

TSFstats<-MetaRareSap %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment, TimeSinceFire) %>%
  summarise(n_obs = n(),
            Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), 
            sd = sd(S.obs), 
            se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstats 


##
##
#******************************************************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS------------------------------------------------------------------------
#******************************************************************************************************************************----
TrtSap <-  MetaRareSap%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>%
  mutate(percent = (mean - first(mean))/first(mean)*100);TrtSap 

TSFburnedSap <- BurnedSap%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>%
  mutate(percent = (mean - first(mean))/first(mean)*100);TSFburnedSap

TSFunburnedSap <- UnburnedSap%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);TSFunburnedSap




# Export files----------------------------------------------------------------------------------------------------------------------
dir.create("Analysis/Diversity/Saprobes/Tables")
dir.create("Analysis/Diversity/Saprobes/Tables/DescriptiveStats")
dir.create("Analysis/Diversity/Saprobes/Tables/PercentChange")

#Descriptive statistics........................................................................................................
capture.output(TrtStats, file="Analysis/Diversity/Saprobes/Tables/DescriptiveStats/Alpha-DesciptiveStats-Trt.csv") 
capture.output(TSFstats, file="Analysis/Diversity/Saprobes/Tables/DescriptiveStats/Alpha-DesciptiveStats-TRT-TSF.csv") 


#Percent change calcualtions ..................................................................................................
capture.output(TrtSap, file="Analysis/Diversity/Saprobes/Tables/PercentChange/Alpha-PercentChange-Trt.csv") 
capture.output(TSFburnedSap, file="Analysis/Diversity/Saprobes/Tables/PercentChange/Alpha-PercentChange-TSF-Burned.csv") 
capture.output(TSFunburnedSap, file="nalysis/Diversity/Saprobes/Tables/PercentChange/Alpha-PercentChange-TSF-Unburned.csv") 



##
##
#******************************************************************************************************************************----
# CHECK DATA FOR NORMALITY---------------------------------------------------------------------------------------------------------
#******************************************************************************************************************************----
#Check data normality to decide between anova and KW test

#ASV data has a normal disribution ................................
par(mfrow=c(3,1))
hist(MetaRareSap$S.obs)
hist(BurnedSap$S.obs)
hist(UnburnedSap$S.obs)

#shapiro test ---Data is not normal, had to normalize
shapiro.test(MetaRareSap$S.obs)#------------
shapiro.test(BurnedSap$S.obs)#--------------
shapiro.test(UnburnedSap$S.obs)#0.359

shapiro.test(MetaRareSap$simpson)#----------
shapiro.test(BurnedSap$simpson)#normal
shapiro.test(UnburnedSap$simpson)#normal


#Transform data to make normal ........................................
Otu<-sqrt(MetaRareSap$S.obs)
OtuB<-log(BurnedSap$S.obs)
OtuUn<-UnburnedSap$S.obs

par(mfrow=c(3,1))
hist(Otu)
hist(OtuB)
hist(OtuUn)

#shapiro test ---Data is normal will use ANOVA to test.................
shapiro.test(Otu)#0.1288
shapiro.test(OtuB)#0.07837
shapiro.test(OtuUn)#0.4308



#
#
#*******************************************************************************************************************************----
# PLOTS SPECIES RICHNESS (ASV's)----------------------------------------------------------------------------------------------------
#*******************************************************************************************************************************----
# TREATMENT -------------------------------------------------------------------------------------------
SppTrtSap<-ggplot(MetaRareSap, aes(x=Treatment, y=S.obs)) + 
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  labs(x = "Treatment", y= " sqrt (Species Richness (ASV's))")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); SppTrtSap

SppTrtSap1<-SppTrtSap + stat_compare_means(method = "anova", label.y = 170, size=6); SppTrtSap1



# SPECIES RICHNESS -- SITES -----------------------------------------------------------------------------
SppSiteSap<-ggplot(MetaRareSap, aes(x=Site, y=S.obs)) + 
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  labs(x = "Site", y= " sqrt (Species Richness (ASV's))")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));SppSiteSap


# SPECIES RICHNESS - TSF  (looks same as old one)-----------------------------------------------------------------------------------------
TrtTSFsap<- ggplot(MetaRareSap, aes(x=TimeSinceFire, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85, alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#45877f","#a2673f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  labs(x = "Time-Since-Fire (days)",
       y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");TrtTSFsap

TrtTSFsap2<-TrtTSFsap+stat_compare_means(method = "anova", label.y=142, label =  "p.signif", 
             show.legend = FALSE, size=6);TrtTSFsap2

# BURNED SPECIES RICHNESS W TIME-SINCE-FIRE --------------------------------------------
TSFsapB<-ggplot(BurnedSap, aes(x=TimeSinceFire, y=S.obs, group=1)) +
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar",
               size=.85, alpha=0.7,position = position_dodge(0.01))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+  
  labs(x = "Time-Since-Fire (days)", 
       y= "Species Richness (ASV's)")+ 
  theme(legend.position="bottom"); TSFsapB

# UNBURNED SPECIES RICHNESS W TIME-SINCE-FIRE ------------------------------------------
TSFsapUn<-ggplot(UnburnedSap, aes(x=TimeSinceFire, y=S.obs, group=1)) +
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", 
               size=.85, alpha=0.7,position = position_dodge(0.01))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+  
  labs(x = "Time-Since-Fire (days)", 
       y= "Species Richness (ASV's)")+ 
  theme(legend.position="bottom"); TSFsapUn





#
#
#**************************************************************************************************************----
#------------------     EXPORT GRAPHS------------------------------------------------------------------------------
#**************************************************************************************************************----
dir.create("Analysis/Diversity/Saprobes/Graphs")

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Trt-TSF.pdf", height=6, width=8)
SppTrtSap1
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Trt-Sites.pdf", height=6, width=8)
SppSiteSap
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Trt-TSF.pdf", height=6, width=8)
TrtTSFsap2
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-TSF-Burned.pdf", height=6, width=8)
TSFsapB
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-TSF-Unburned.pdf", height=6, width=8)
TSFsapUn
dev.off()

#
###
#####
#**********************************************************************************************************----
# ------------------------------      NESTED MODEL SELECTION  -------------------------------------------------
#**********************************************************************************************************----
# TREATMENT .......................................................
# *
Otu# Metarare
OtuB#Burned
OtuUn#Unburned

attach(MetaRareSap)
library(lmerTest)


#First check the level of nestedness to see which one to use ..........................................
# Looks like OTUlm1 is the model to use 
OTUlm1<- lmer(Otu ~ Treatment + (1|Site/Plot/Subplot), data = MetaRareSap, REML=TRUE);summary(OTUlm1)
OTUlm2<- lmer(Otu ~ Treatment + (1|Site/Plot), data = MetaRareSap, REML=TRUE);summary(OTUlm2)#0.0683
OTUlm3<- lmer(Otu ~ Treatment + (1|Site), data = MetaRareSap, REML=TRUE);summary(OTUlm3)#0.000553
OTUlm4<- lmer(Otu ~ Treatment + (1|Plot), data = MetaRareSap, REML=TRUE);summary(OTUlm4)#0.0762

anova(OTUlm1,OTUlm2)#Mod2 better
anova(OTUlm2,OTUlm3)#Mod2 better
anova(OTUlm2,OTUlm4)#Mod4 better

anova(OTUlm1,OTUlm2, OTUlm3,OTUlm4)#Mod 2 better


#
###
#####
#******************************************************************************************************************----
# ------------------------------      TEST VARIABLE SIGNIFICANCE  ---------------------------------------------
#******************************************************************************************************************----

#Effect of treatment on species richness (#Uses transform data......................................................
OtuTrtTsfLM<- lmer(Otu ~ Treatment * TimeSinceFire + (1|Plot), data = MetaRareSap, REML=TRUE);summary(OtuTrtTsfLM)
OtuTrtTsfAV<-anova(OtuTrtTsfLM);OtuTrtTsfAV#all significant 


#Final nestedness test .................................................. .........................................
OtuTrtLM<- lmer(Otu ~ Treatment + (1|Plot), data = MetaRareSap, REML=TRUE);summary(OtuTrtLM)#0.0762 not sign
OtuTrtAV<-anova(OtuTrtLM);OtuTrtAV#Gives us the same pvalue 

detach(MetaRareSap)
#Test TimeSinceFire but has to be done independently for burned and unburned sites

#
##
#*********************************************************************************************************************************----
# ** UNBURNED ONLY -------------------------------------------------------------------------------------------------------------------
#*********************************************************************************************************************************----

#TIME SINCE FIRE 
attach(Burned)
OTUlm1B<- lmer(OtuB ~ TimeSinceFire + (1|Site/Plot/Subplot), data = BurnedSap, REML=TRUE);summary(OTUlm1B)#does not work
OTUlm2B<- lmer(OtuB ~ TimeSinceFire + (1|Site/Plot), data = BurnedSap, REML=TRUE);summary(OTUlm2B)#no converg
OTUlm3B<- lmer(OtuB ~ TimeSinceFire + (1|Site), data = BurnedSap, REML=TRUE);summary(OTUlm3B)#
OTUlm4B<- lmer(OtuB ~ TimeSinceFire + (1|Plot), data = BurnedSap, REML=TRUE);summary(OTUlm4B)#

anova(OTUlm2B,OTUlm3B)#Mod3 better
anova(OTUlm3B,OTUlm4B)#Same AIC
anova(OTUlm2B, OTUlm3B,OTUlm4B)#Mod 3 and 4 are the same, will proceed with plot to keep analysis same as above


#TEST FOR SIGNIFICANCE ...................................................................................................
OtuBtsfLM<- lmer(OtuB ~ TimeSinceFire + (1|Plot), data = BurnedSap, REML=TRUE);summary(OtuBtsfLM)#Intercept only one sign
OtuBtsfAV<-anova(OtuBtsfLM);OtuBtsfAV#Pvalue = 0.6704

#Post-hoc analysis not needed since the reg is not significant 
detach(Burned)



#*********************************************************************************************************************************----
# ** UNBURNED ONLY -------------------------------------------------------------------------------------------------------------------
#*********************************************************************************************************************************----
#TIME SINCE FIRE 
attach(Unburned)
OTUlm1Un<- lmer(OtuUn ~ TimeSinceFire + (1|Site/Plot/Subplot), data = UnburnedSap, REML=TRUE);summary(OTUlm1Un)#does not work
OTUlm2Un<- lmer(OtuUn ~ TimeSinceFire + (1|Site/Plot), data = UnburnedSap, REML=TRUE);summary(OTUlm2Un)#no converg
OTUlm3Un<- lmer(OtuUn ~ TimeSinceFire + (1|Site), data = UnburnedSap, REML=TRUE);summary(OTUlm3Un)#
OTUlm4Un<- lmer(OtuUn ~ TimeSinceFire + (1|Plot), data = UnburnedSap, REML=TRUE);summary(OTUlm4Un)#

anova(OTUlm2Un,OTUlm3Un)#Mod2 better
anova(OTUlm2Un,OTUlm4Un)#Mod 4 better, proceed w 4
anova(OTUlm2Un, OTUlm3Un,OTUlm4Un)#Mod 4


#TEST FOR SIGNIFICANCE .....................................................................................................
OtuUntsfLM<- lmer(OtuUn ~ TimeSinceFire + (1|Plot), data = UnburnedSap, REML=TRUE);summary(OtuUntsfLM)#Intercept & 11 sign
OtuUntsfAV<-anova(OtuUntsfLM);OtuUntsfAV#Pvalue = 0.007247

#Post-hoc analysis not needed since the reg is not significant 
detach(Unburned)



##
##
#********************************************************************************************************************************----
#EXPORT RESULTS FROM ANOVA-----------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
#Treatment * TSF
capture.output(summary(OtuTrtTsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-TSF-Spp.CSV")
capture.output(OtuTrtTsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-TSF-Spp.CSV")

#Interaction of tretment 
capture.output(summary(OtuTrtLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp.CSV")
capture.output(OtuTrtAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp.CSV")

#BurnedSap
capture.output(summary(OtuBtsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp-Burned.CSV")
capture.output(OtuBtsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp-Burned.CSV")

##UnburnedSap
capture.output(summary(OtuUntsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp-Unburned.CSV")
capture.output(OtuUntsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp-Unburned.CSV")






#
#
#*****************************************************************************************************************----
#--------------------    SAME CALCULATIONS AS ABOVE USING SIMPSONS DIVERSITY INDEX -----------------------------------
#*****************************************************************************************************************----

# TREATMENT -------------------------------------------------------------------------------------------
attach(MetaRareSap)

#Not normal
shapiro.test(MetaRareSap$simpson)# 
shapiro.test(BurnedSap$simpson)#
shapiro.test(UnburnedSap$simpson)#

shapiro.test(log2(MetaRareSap$simpson))# 
shapiro.test(log2(BurnedSap$simpson))#
shapiro.test(log2(UnburnedSap$simpson))#



SimpSapTrt<-ggplot(MetaRareSap, aes(x=Treatment, y=simpson)) +   
              geom_boxplot(aes(fill=Treatment))+ 
              theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+
              labs(x = "", y= "Simpson's Index")+ theme(legend.position="bottom")+
              scale_fill_manual(values=c("#45877f","#a2673f")); SimpSapTrt

SimpSapTrt1<-SimpSapTrt + stat_compare_means(method = "kruskal", label.y = 1);SimpSapTrt1

# SPECIES RICHNESS -- SITES -----------------------------------------------------------------------------
SimpSapSites<-ggplot(MetaRareSap, aes(x=Site, y=simpson)) +   
                geom_boxplot(aes(fill=Treatment))+ 
                theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+   
              labs(x = "", y= "Simpson's Index")+theme(legend.position="bottom")+
              scale_fill_manual(values=c("#45877f","#a2673f"));SimpSapSites


# SPECIES RICHNESS - TIMES SINCE FIRE -----------------------------------------------------------------------------------------
SimpSamTrtTSF<- ggplot(MetaRareSap, aes(x=TimeSinceFire, y=simpson, group=Treatment, col=Treatment))+
                  stat_summary(fun=mean,geom="line", size=1.3)+
                  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
                         alpha=0.7,position = position_dodge(0.01))+
                  scale_color_manual(values=c("#45877f","#a2673f")) + 
                  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", fill=NA, size=.7),
                         text = element_text(size=23),
                         axis.text.y = element_text(size=19,hjust=0.5),
                         axis.text.x = element_text(size=19))+
                  labs(x = "Time Since Fire (Years)", y= "Simpson's Index")+
                  theme(legend.position="bottom");SimpSamTrtTSF

SimpSapTrtTSF2<-SimpSamTrtTSF +stat_compare_means(method = "kruskal", label.y=0.965,
                    label =  "p.signif",how.legend = FALSE);SimpSapTrtTSF2

# BURNED SPECIES RICHNESS W TIME-SINCE-FIRE --------------------------------------------
TSFsapBsimp<-ggplot(BurnedSap, aes(x=TimeSinceFire, y=simpson, group=1)) +
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", 
               size=.85, alpha=0.7,position = position_dodge(0.01))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+  
  labs(x = "Time-Since-Fire (days)", 
       y= "Species Richness (ASV's)")+ 
  theme(legend.position="bottom"); TSFsapBsimp

# UNBURNED SPECIES RICHNESS W TIME-SINCE-FIRE ------------------------------------------
TSFsapUnSimp<-ggplot(UnburnedSap, aes(x=TimeSinceFire, y=simpson, group=1)) +
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", 
               size=.85, alpha=0.7,position = position_dodge(0.01))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+  
  labs(x = "Time-Since-Fire (days)", 
       y= "Species Richness (ASV's)")+ 
  theme(legend.position="bottom"); TSFsapUnSimp




#
#**************************************************************************************************************----
#EXPORT PLOTS -----------------------------------------------------------------------------------------------
#**************************************************************************************************************----
dir.create("Analysis/Diversity/Saprobes/Graphs")

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Treatment.pdf", height=6, width=8)
SimpSapTrt1
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Trt-Sites.pdf", height=6, width=8)
SimpSapSites
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/SPP-Trt-TSF.pdf", height=6, width=8)
SimpSapTrtTSF2
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/Simpson-SPP-TSF-Burned.pdf", height=6, width=8)
TSFsapBsimp
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/Simpson-SPP-TSF-Unburned.pdf", height=6, width=8)
TSFsapUnSimp
dev.off()





#
###
#####
#**********************************************************************************************************----
# ------------------------------      NESTED MODEL SELECTION  -------------------------------------------------
#**********************************************************************************************************----
# TREATMENT .......................................................
# *

attach(MetaRareSap)
library(lmerTest)


#First check the level of nestedness to see which one to use ..........................................
# Looks like simpsonlm1 is the model to use 
simpsonlm1<- lmer(simpson ~ Treatment + (1|Site/Plot/Subplot), data = MetaRareSap, REML=TRUE);summary(simpsonlm1)
simpsonlm2<- lmer(simpson ~ Treatment + (1|Site/Plot), data = MetaRareSap, REML=TRUE);summary(simpsonlm2)#0.0683
simpsonlm3<- lmer(simpson ~ Treatment + (1|Site), data = MetaRareSap, REML=TRUE);summary(simpsonlm3)#0.000553
simpsonlm4<- lmer(simpson ~ Treatment + (1|Plot), data = MetaRareSap, REML=TRUE);summary(simpsonlm4)#0.0762

anova(simpsonlm1,simpsonlm2)#Mod2 better
anova(simpsonlm2,simpsonlm3)#Mod2 better
anova(simpsonlm2,simpsonlm4)#Mod4 better
anova(simpsonlm1,simpsonlm2, simpsonlm3,simpsonlm4)#Mod 4 better


#
###
#####
#******************************************************************************************************************----
# ------------------------------      TEST VARIABLE SIGNIFICANCE  ---------------------------------------------
#******************************************************************************************************************----

#Effect of treatment on species richness (#Uses transform data......................................................
simpsonTrtTsfLM<- lmer(simpson ~ Treatment * TimeSinceFire + (1|Plot), data = MetaRareSap, REML=TRUE);summary(simpsonTrtTsfLM)
simpsonTrtTsfAV<-anova(simpsonTrtTsfLM);simpsonTrtTsfAV#Treatment not significant, only TSF

#Final nestedness test .................................................. .........................................
simpsonTrtLM<- lmer(simpson ~ Treatment + (1|Plot), data = MetaRareSap, REML=TRUE);summary(simpsonTrtLM)#0.0762 not sign
simpsonTrtAV<-anova(simpsonTrtLM);simpsonTrtAV#Gives us the same pvalue 

detach(MetaRareSap)
#Test TimeSinceFire but has to be done independently for burned and unburned sites

#
##
#*********************************************************************************************************************************----
# ** UNBURNED ONLY -------------------------------------------------------------------------------------------------------------------
#*********************************************************************************************************************************----

#TIME SINCE FIRE 
attach(Burned)
simpsonlm1B<- lmer(simpson ~ TimeSinceFire + (1|Site/Plot/Subplot), data = BurnedSap, REML=TRUE);summary(simpsonlm1B)#does not work
simpsonlm2B<- lmer(simpson ~ TimeSinceFire + (1|Site/Plot), data = BurnedSap, REML=TRUE);summary(simpsonlm2B)#no converg
simpsonlm3B<- lmer(simpson ~ TimeSinceFire + (1|Site), data = BurnedSap, REML=TRUE);summary(simpsonlm3B)#
simpsonlm4B<- lmer(simpson ~ TimeSinceFire + (1|Plot), data = BurnedSap, REML=TRUE);summary(simpsonlm4B)#mo convergence

anova(simpsonlm2B,simpsonlm3B)#Mod2 better
anova(simpsonlm2B,simpsonlm4B)#Mod 2 better
anova(simpsonlm2B, simpsonlm3B,simpsonlm4B)#M as above
#Mod 3 is the only one that converges will proceed with site as model to use 

#TEST FOR SIGNIFICANCE ...................................................................................................
simpsonBtsfLM<- lmer(simpson~ TimeSinceFire + (1|Site), data = BurnedSap, REML=TRUE);summary(simpsonBtsfLM)#Intercept only one sign
simpsonBtsfAV<-anova(simpsonBtsfLM);simpsonBtsfAV#Pvalue = 0.7918

#Post-hoc analysis not needed since the reg is not significant 
detach(Burned)



#*********************************************************************************************************************************----
# ** UNBURNED ONLY -------------------------------------------------------------------------------------------------------------------
#*********************************************************************************************************************************----
#TIME SINCE FIRE 
attach(Unburned)
simpsonlm1Un<- lmer(simpson ~ TimeSinceFire + (1|Site/Plot/Subplot), data = UnburnedSap, REML=TRUE);summary(simpsonlm1Un)#does not work
simpsonlm2Un<- lmer(simpson ~ TimeSinceFire + (1|Site/Plot), data = UnburnedSap, REML=TRUE);summary(simpsonlm2Un)#no converg
simpsonlm3Un<- lmer(simpson ~ TimeSinceFire + (1|Site), data = UnburnedSap, REML=TRUE);summary(simpsonlm3Un)#
simpsonlm4Un<- lmer(simpson ~ TimeSinceFire + (1|Plot), data = UnburnedSap, REML=TRUE);summary(simpsonlm4Un)#no convergence

anova(simpsonlm2Un,simpsonlm3Un)#Mod2 better
anova(simpsonlm2Un,simpsonlm4Un)#Mod 4 better, proceed w 4
anova(simpsonlm2Un, simpsonlm3Un,simpsonlm4Un)#Mod 4


#TEST FOR SIGNIFICANCE .....................................................................................................
simpsonUntsfLM<- lmer(simpsonUn ~ TimeSinceFire + (1|Plot), data = UnburnedSap, REML=TRUE);summary(simpsonUntsfLM)#Intercept & 11 sign
simpsonUntsfAV<-anova(simpsonUntsfLM);simpsonUntsfAV#Pvalue = 0.007247

#Post-hoc analysis not needed since the reg is not significant 
detach(Unburned)



##
##
#********************************************************************************************************************************----
#EXPORT RESULTS FROM ANOVA-----------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
#Treatment * TSF
capture.output(summary(simpsonTrtTsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-TSF-Spp.CSV")
capture.output(simpsonTrtTsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-TSF-Spp.CSV")

#Interaction of tretment 
capture.output(summary(simpsonTrtLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp.CSV")
capture.output(simpsonTrtAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp.CSV")

#BurnedSap
capture.output(summary(simpsonBtsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp-Burned.CSV")
capture.output(simpsonBtsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp-Burned.CSV")

##UnburnedSap
capture.output(summary(simpsonUntsfLM), file="Analysis/Diversity/Saprobes/Tables/LM-Trt-Spp-Unburned.CSV")
capture.output(simpsonUntsfAV, file = "Analysis/Diversity/Saprobes/Tables/AnovaLM-Trt-Spp-Unburned.CSV")







































































