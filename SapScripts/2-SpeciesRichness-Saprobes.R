#Date last ran June 25, 2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY--------------------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")
library(ggplot2)
library(ggpubr)
library(dplyr)


#Load Data--------------------------------------------------------------------------------------------------------
MetaRareSap<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.csv", na.strings = "N/A", header = TRUE) 
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
MetaRareSap$Treatment <- try(relevel(MetaRareSap$Treatment , "Unburned"));levels(MetaRareSap$Treatment)

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
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TrtStats # B=8.413793, U=15.514286

TSFstats<-MetaRareSap %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment, TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstats 


##
##
#******************************************************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS------------------------------------------------------------------------
#******************************************************************************************************************************----
TrtSap <-  MetaRareSap%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);TrtSap 

TSFburnedSap <- BurnedSap%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
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
shapiro.test(MetaRareSap$S.obs)
shapiro.test(BurnedSap$S.obs)#0.0303
shapiro.test(UnburnedSap$S.obs)#0.359

shapiro.test(MetaRareSap$simpson)
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

#shapiro test ---Data is normal will use ANOVA to test
shapiro.test(Otu)#0.1288
shapiro.test(OtuB)#0.07837
shapiro.test(OtuUn)#0.4308



#
#
#*******************************************************************************************************************************----
# PLOTS SPECIES RICHNESS (ASV's)----------------------------------------------------------------------------------------------------
#*******************************************************************************************************************************----


# TREATMENT -------------------------------------------------------------------------------------------
SppTrtSap<-ggplot(MetaRareSap, aes(x=Treatment, y=Otu)) +   
            geom_boxplot(aes(fill=Treatment))+ 
            theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=.7), text = element_text(size=23),
              axis.text.y = element_text(size=19,hjust=0.5),axis.text.x = element_text(size=19))+   
              labs(x = "Treatment", y= " sqrt (Species Richness (ASV's))")+theme(legend.position="bottom")+
              scale_fill_manual(values=c("#a2673f", "#45877f")); SppTrtSap

SppTrtSap1<-SppTrtSap + stat_compare_means(method = "anova", label.y = 14,); SppTrtSap1



# SPECIES RICHNESS -- SITES -----------------------------------------------------------------------------
SppSiteSap<-ggplot(MetaRareSap, aes(x=Site, y=Otu)) +   
                geom_boxplot(aes(fill=Treatment))+ 
                theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),axis.text.x = element_text(size=19))+   
                labs(x = "Sites", y=" sqrt (Species Richness (ASV's))")+theme(legend.position="bottom")+
                scale_fill_manual(values=c("#a2673f", "#45877f"));SppSiteSap


# SPECIES RICHNESS - TSF  (looks same as old one)-----------------------------------------------------------------------------------------
TrtTSFsap<- ggplot(MetaRareSap, aes(x=TimeSinceFire, y=Otu, group=Treatment, col=Treatment))+
              stat_summary(fun=mean,geom="line", size=1.3)+
              stat_summary(fun.data = mean_se,geom = "errorbar", size=.85, alpha=0.7,
                           position = position_dodge(0.01))+
              scale_color_manual(values=c("#a2673f", "#45877f")) + 
              theme_bw()+theme(panel.grid.major = element_blank()
                               ,panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                    text = element_text(size=23),
                    axis.text.y = element_text(size=19,hjust=0.5),
                    axis.text.x = element_text(size=19))+   
              labs(x = "Time-Since-Fire (days)", y= " sqrt ((Species Richness (ASV's))")+
              theme(legend.position="bottom");TrtTSFsap

TrtTSFsap2<-TrtTSFsap+stat_compare_means(method = "anova", label.y=12, label =  "p.signif", 
                    show.legend = FALSE, size=10);TrtTSFsap2

# BURNED SPECIES RICHNESS W TIME-SINCE-FIRE --------------------------------------------
TSFsapB<-ggplot(BurnedSap, aes(x=TimeSinceFire, y=S.obs, group=1)) +
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
#
#***********************************************************************************************************************----
#TEST SIGNIFICANCE --------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
# TREATMENT -------------Kruskal Wallis-treatment-all data (no subsetting)
attach(MetaRareSap)
SapTrtLm <- lm(Otu ~ Treatment, data = MetaRareSap);SapTrtLm
AovTrt<- aov(SapTrtLm);
AovSumTrt<-summary(AovTrt);AovSumTrt#0.000566 ***
TukTrt<- TukeyHSD(AovTrt);TukTrt
detach(MetaRareSap)

#TIME-SINCE-FIRE UnburnedSap DATA ------------------------------------------------------------------
attach(UnburnedSap)
SapTsfUnLm <- lm(OtuUn ~ TimeSinceFire, data = UnburnedSap);SapTsfUnLm 
AovTsfUn<- aov(SapTsfUnLm );
AovSumTSFUn<-summary(AovTsfUn);AovSumTSFUn#1.94e-12 ***
TukTsfUn<- TukeyHSD(AovTsfUn);TukTsfUn#11-2, 11-3, 11-5
detach(UnburnedSap)

#TIME-SINCE-FIRE BurnedSap DATA --------------------------------------------------------------------
attach(BurnedSap)
SapTsfBLm <- lm(OtuB ~ TimeSinceFire, data = BurnedSap);SapTsfBLm 
AovTsfB<- aov(SapTsfBLm );
AovSumTSFB<-summary(AovTsfB);AovSumTSFB#0.435------------------------NOT SIGN
TukTsfB<- TukeyHSD(AovTsfB);TukTsfB##No SIgn
detach(BurnedSap)


##
##
#**************************************************************************************************************----
#EXPORT RESULTS FROM ANOVA--------------------------------------------------------------------------------
#**************************************************************************************************************----
capture.output(AovSumTrt, file="Analysis/Diversity/Saprobes/Tables/TreatmentAnova.CSV")
capture.output(TukTrt, file = "Analysis/Diversity/Saprobes/Tables/Treatment-TukeyTest.CSV")

#UnburnedSap
capture.output(AovSumTSFUn, file="Analysis/Diversity/Saprobes/Tables/TSF-Anova-UnburnedSap.CSV")
capture.output(TukTsfUn, file = "Analysis/Diversity/Saprobes/Tables/TSF-TukeyTest-UnburnedSap.CSV")

#BurnedSap
capture.output(AovSumTSFB, file="Analysis/Diversity/Saprobes/Tables/TSF-Anova-BurnedSap.CSV")
capture.output(TukTsfB, file = "Analysis/Diversity/Saprobes/Tables/TSF-TukeyTest-BurnedSap.CSV")





#
#
#*****************************************************************************************************************----
#--------------------    SAME CALCULATIONS AS ABOVE USING SIMPSONS DIVERSITY INDEX -----------------------------------
#*****************************************************************************************************************----

# TREATMENT -------------------------------------------------------------------------------------------
attach(MetaRareSap)

shapiro.test(MetaRareSap$simpson)# cannot transform

shapiro.test(BurnedSap$simpson)# normal
shapiro.test(UnburnedSap$simpson)# normal


SimpSapTrt<-ggplot(MetaRareSap, aes(x=Treatment, y=simpson)) +   
              geom_boxplot(aes(fill=Treatment))+ 
              theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+
              labs(x = "", y= "Simpson's Index")+ theme(legend.position="bottom")+
              scale_fill_manual(values=c("#45877f","#a2673f")); SimpSapTrt

SimpSapTrt1<-SimpSapTrt + stat_compare_means(method = "kruskal", label.y = 50,);SimpSapTrt1

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

SimpSapTrtTSF2<-SimpSamTrtTSF +stat_compare_means(method = "kruskal", label.y=30,
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
#
#*************************************************************************************************************----
#TEST SIGNIFICANCE -----------------------------------------------------------------------------------------------
#*************************************************************************************************************----

# TREATMENT -------------
attach(MetaRareSap)
SapTrtLmS <- lm(simpson ~ Treatment, data = MetaRareSap);SapTrtLmS
AovTrtS<- aov(SapTrtLmS);
AovSumTrtS<-summary(AovTrtS);AovSumTrtS#0.000519 ***
detach(MetaRareSap)

#TIM-SINCE-FIRE UnburnedSap DATA ------------------------------------------------------------------------------------
attach(UnburnedSap)
SapTsfUnLmS <- lm(simpson ~ TimeSinceFire, data = UnburnedSap);SapTsfUnLmS 
AovTsfUnS <- aov(SapTsfUnLmS);
AovSumTSFUnS <- summary(AovTsfUnS);AovSumTSFUnS#0.000161*
TukTsfUnS <- TukeyHSD(AovTsfUnS);TukTsfUnS#11-2, 5-2, 11-3
detach(UnburnedSap)

#TIME-SINCE-FIRE BurnedSap DATA -------------------------------------------------------------------------------------
attach(BurnedSap)
SapTsfBLmS <- lm(simpson ~ TimeSinceFire, data = BurnedSap);SapTsfBLmS 
AovTsfBS<- aov(SapTsfBLmS);AovTsfBS
AovSumTSFBS<-summary(AovTsfBS);AovSumTSFBS#0.333 -----NOT Significant
TukTsfBS <- TukeyHSD(AovTsfBS);TukTsfBS # no sign
detach(BurnedSap)



#**************************************************************************************************************----
#EXPORT RESULTS FROM ANOVA--------------------------------------------------------------------------------
#**************************************************************************************************************----
capture.output(AovSumTrtS, file="Analysis/Diversity/Saprobes/Tables/Simpson-TreatmentAnova.CSV")
capture.output(TukTrtS, file = "Analysis/Diversity/Saprobes/Tables/Simpson-Treatment-TukeyTest.CSV")

#UnburnedSap
capture.output(AovSumTSFUnS, file="Analysis/Diversity/Saprobes/Tables/Simpson-TSF-Anova-UnburnedSap.CSV")
capture.output(TukTsfUnS, file = "Analysis/Diversity/Saprobes/Tables/Simpson-TSF-TukeyTest-UnburnedSap.CSV")

#BurnedSap
capture.output(AovSumTSFBS, file="Analysis/Diversity/Saprobes/Tables/SImpson-TSF-Anova-BurnedSap.CSV")
capture.output(TukTsfBS, file = "Analysis/Diversity/Saprobes/Tables/Simpson-TSF-TukeyTest-BurnedSap.CSV")

























































