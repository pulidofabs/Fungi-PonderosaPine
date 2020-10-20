#10-9-2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")


#Load libraries---------
library(ggplot2)
library(ggpubr)
library(data.table)# to calculate percent change
library(dplyr)


#Load Data--------------------------------------------------------------------------
MetaRare<- read.csv("Analysis/Metadata/Fungi/MetaRareCoreFun.csv", na.strings = "N/A", header = TRUE) #Rarefied metadata


# Convert to factor------------------------------------------------
MetaRare$FireYear <- as.factor(MetaRare$FireYear)
MetaRare$TimeSinceFire <- as.factor(MetaRare$TimeSinceFire)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRare)
Unburned<-MetaRare[which(Treatment== "Unburned"), ]
head(Unburned)
dim(Unburned)#106x33

# * * BURNED PLOTS ------------------------------
Burned<-MetaRare[which(Treatment== "Burned"), ]
head(Burned)
dim(Burned)#104x33


#
#
#***********************************************************************************************----
# ALPHA DIVERSITY- MEAN CALCULATIONS----------------------------------------------------------------
#***********************************************************************************************----

MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment <- try(relevel(MetaRare$Treatment , "Unburned"))


TrtStats<-MetaRare %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TrtStats # B=203.6827, U=232.7925

TSFstats<-MetaRare %>%
  filter(!is.na(S.obs)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstats # 2-205-2963, 3=204.7925, 5= 216.6604, 11=248.7200 


TSFstatsB<-Burned %>%
  filter(!is.na(S.obs)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstatsB # 2-202.5185, 3=228.6538, 5= 199.6923, 11=183.1200


TSFstatsUn<-Unburned %>%
  filter(!is.na(S.obs)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstatsUn # 2-208.0741, 3=181.8148, 5= 233.0000, 11=314.3200


#eXPORT Restults........................................................................................................
dir.create("Analysis/Diversity/Fungi/Tables/")
dir.create("Analysis/Diversity/Fungi/Tables/DescriptiveStats")

write.csv(TrtStats, "Analysis/Diversity/Fungi/Tables/DescriptiveStats/TRT-DesStats-Fungi.csv")
write.csv(TSFstats, file = "Analysis/Diversity/Fungi/Tables/DescriptiveStats/TSF-DesStats-Fungi.csv")
write.csv(TSFstatsB, file = "Analysis/Diversity/Fungi/Tables/DescriptiveStats/TSF-DesStats-Fungi-Burned.csv")
write.csv(TSFstatsUn, file = "Analysis/Diversity/Fungi/Tables/DescriptiveStats/TSF-DesStats-Fungi-Unburned.csv")


#
#
#***************************************************************************************************************----
# TRETMENT ALPHA PERCENT CHANGE-------------------------------------------------------------------------------------
#***************************************************************************************************************----
FungiTrt <-  MetaRare%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);FungiTrt #-12.5

FungiTsfBurned <- Burned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);FungiTsfBurned

FungiTsfUnburned <- Unburned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);FungiTsfUnburned



#EXPORTFILES---------------------------------------------------------------------------------------------------------
dir.create("Analysis/Diversity/Fungi/Tables/PerChange")

write.csv(FungiTrt, "Analysis/Diversity/Fungi/Tables/PerChange/TRT-PerChange-Fungi.csv")
write.csv(FungiTsfBurned, file = "Analysis/Diversity/Fungi/Tables/PerChange/TSF-PerChange-Fungi-Burned.csv")
write.csv(FungiTsfUnburned, file = "Analysis/Diversity/Fungi/Tables/PerChange/TSF-PerChange-Fungi-Unburned.csv")




#*****************************************************************************************************************----
# ---------------------- QUALITY CHECK -------------------------------------------------------------------------------
#*****************************************************************************************************************----
# CHECK DATA FOR NORMALITY--------------------------------------------------------------------------
par(mfrow=c(3,1))
hist(MetaRare$S.obs)
hist(Burned$S.obs)
hist(Unburned$S.obs)

#shapiro test 
shapiro.test(MetaRare$S.obs)#0.0002785
shapiro.test(Burned$S.obs)#0.0163
shapiro.test(Unburned$S.obs)#0.04436


#GLM models dependent variable Yi does NOT need to be normally distributed
#lmer models require data to be normally distributed


#tranform to see if it makes a difference, if yes = ANOVA if no WILCOXON TEST
S.obs<-log(MetaRare$S.obs)
S.obsB<-log(Burned$S.obs)
S.obsU<-log(Unburned$S.obs)

#Log transformed normalized data
shapiro.test(S.obs)#0.6245
shapiro.test(S.obsB)#0.6721
shapiro.test(S.obsU)#0.1395


#
#
#-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*_*_*_*
# PLOTS SPECIES RICHNESS (ESV's)---------------------------------------------------------------------
#-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*_*_*_*

# TREATMENT -------------------------------------------------------------------------------------------
SppTrt<-ggplot(MetaRare, aes(x=Treatment, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Sites", y= "Species Richness (ASV's)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); SppTrt

SppTrt1<-SppTrt + stat_compare_means(method = "anova", size=6, 
              label.x = 2, label.y = 400 ); SppTrt1

# SPECIES RICHNESS -- SITES --------------------------------------------------------------------------------------
SppTSFBox<-ggplot(MetaRare, aes(x=Site, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  labs(x = "Sites", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));SppTSFBox



# SPECIES RICHNESS - TSF -----------------------------------------------------------------------------------------
SppTSF<-ggplot(MetaRare, aes(x=TimeSinceFire, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
               alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#45877f","#a2673f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Time Since Fire (Years)", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");SppTSF

SppTSF1<-SppTSF+stat_compare_means(method = "anova", label.y=350, label =  "p.signif", 
                  show.legend = FALSE, size = 6);SppTSF1



SppFY<-ggplot(MetaRare, aes(x=FireYear, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
               alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#45877f","#a2673f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Time Since Fire (Years)", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");SppFY

SppFY1<-SppFY+stat_compare_means(method = "anova", label.y=350, label =  "p.signif", 
                show.legend = FALSE, size = 6);SppFY1






# EXPORT PLOTS ---------------------------------------------------------------------------------------------------------
dir.create(file.path("Analysis/Diversity/Fungi/Graphs/"), recursive = TRUE)

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-Treatment.pdf", height=6, width=8)
SppTrt1
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-TSF-Boxplot.pdf", height=6, width=8)
SppTSFBox
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-TSF-Line.pdf", height=8, width=10)
SppTSF1
dev.off()


pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-FY-Line.pdf", height=8, width=10)
SppFY1
dev.off()




#*******************************************************************************************************----
#---------CREATE A PANEL----------------------------------------------------------------------------------
#*******************************************************************************************************----

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-panels.pdf", height=8, width=11, onefile = FALSE)
ggarrange(SppTrt1, SppTSFBox,  ncol=2, nrow=1, labels= c("a","b"), 
          font.label = list(size = 16,face = "bold", color ="black"),
          legend = "bottom",common.legend = TRUE,  align = "hv")
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-panels-TRT-TSF.pdf", height=8, width=11, onefile = FALSE)
ggarrange(SppTrt1, SppTSF1,  ncol=2, nrow=1, labels= c("a","b"), 
          font.label = list(size = 16,face = "bold", color ="black"),
          legend = "bottom",common.legend = TRUE,  align = "hv")
dev.off()





#
#
#*************************************************************************************************************----
#TEST SIGNIFICANCE -----------------------------------------------------------------------------------------------
#*************************************************************************************************************----
# TREATMENT -------------Kruskal Wallis-treatment-all data (no subsetting)
# *** if use Site/Plot model fails to converge
# *** If use only Site model works
# 
attach(MetaRare)
library(lmerTest)

S.obs<-log(MetaRare$S.obs)
S.obsB<-log(Burned$S.obs)
S.obsU<-log(Unburned$S.obs)


#First check the level of nestedness to see which one to use .............................................
# Looks like OTUlm1 is the model to use 
library(lme4)
library(lmerTest)

OTUlm1<- lmer(S.obs ~ Treatment + (1|Site/Plot/Subplot), data = MetaRare, REML=TRUE);summary(OTUlm1)#does not converge
OTUlm2<- lmer(S.obs ~ Treatment + (1|Site/Plot), data = MetaRare, REML=TRUE);summary(OTUlm2)#Model does not converge
OTUlm3<- lmer(S.obs  ~ Treatment + (1|Site), data = MetaRare, REML=TRUE);summary(OTUlm3)#...........best model
OTUlm4<- lmer(S.obs  ~ Treatment + (1|Plot), data = MetaRare, REML=TRUE);summary(OTUlm4)
anova(OTUlm1,OTUlm2)#Mod2 better
anova(OTUlm2,OTUlm3)#Mod2 better
anova(OTUlm2,OTUlm4)#Mod4 better

#Model gave the same results so will proceed with Site level nestedness for all analysis


##
####
#*************************************************************************************************************----
#-------------------  TEST FOR SIGNIFICANCE OF THE DATA ----------------------------------------------------------
#*************************************************************************************************************----
#Effect of treatment on species richness (using transform species richness for significance).........................
OTUlmTrtTSF1<- lmer(S.obs ~ Treatment * TimeSinceFire + (1|Plot), data = MetaRare, REML=TRUE);summary(OTUlmTrtTSF1)
OTUavTrtTSF<-anova(OTUlmTrtTSF1);OTUavTrtTSF##All significant


#Final nestedness test .................................................. ..........................
OTUlm1<- lmer(S.obs~ Treatment + (1|Plot), data = MetaRare, REML=TRUE);summary(OTUlm1)#not sign
OTUav<-anova(OTUlm1);OTUav#not sign 0.1012

#Time since fire will be tested independently between burned and unburned sites







#TIME-SINCE-FIRE BURNED DATA ---------------------------------------------------------------------------------------------------
attach(Burned)
OTUlm1<- lmer(S.obsB ~ TimeSinceFire + (1|Site/Plot/Subplot), data = Burned, REML=TRUE);summary(OTUlm1)#does not converge
OTUlm2<- lmer(S.obsB ~ TimeSinceFire + (1|Site/Plot), data = Burned, REML=TRUE);summary(OTUlm2)#Model does not converge
OTUlm3<- lmer(S.obsB  ~ TimeSinceFire + (1|Site), data = Burned, REML=TRUE);summary(OTUlm3)#Model does not converge
OTUlm4<- lmer(S.obsB  ~ TimeSinceFire + (1|Plot), data = Burned, REML=TRUE);summary(OTUlm4)

anova(OTUlm2,OTUlm3)#Mod3 better but it did not converge
anova(OTUlm3,OTUlm4)#Mod4 better......................................onlyone that converge

#TSF not significant
BurnOtuLm1<- lmer(S.obsB  ~ TimeSinceFire + (1|Plot), data = Burned, REML=TRUE);summary(BurnOtuLm1)#only intercept significant
BurnOTUav<-anova(BurnOtuLm1);BurnOTUav#Pvalue = 0.2075 

#Post-hoc analysis...................................................
#** Not needed since the reg is not significant 
detach(Burned)



#TIM-SINCE-FIRE UNBURNED DATA ------------------------------------------------------------------------------------
attach(Unburned)
OTUlm1<- lmer(S.obsU ~ TimeSinceFire + (1|Site/Plot/Subplot), data = Unburned, REML=TRUE);summary(OTUlm1)#does not converge
OTUlm2<- lmer(S.obsU ~ TimeSinceFire + (1|Site/Plot), data = Unburned, REML=TRUE);summary(OTUlm2)#Model does not converge
OTUlm3<- lmer(S.obsU  ~ TimeSinceFire + (1|Site), data = Unburned, REML=TRUE);summary(OTUlm3)#Model does not converge
OTUlm4<- lmer(S.obsU  ~ TimeSinceFire + (1|Plot), data = Unburned, REML=TRUE);summary(OTUlm4)

anova(OTUlm2,OTUlm3)#Mod3 better but it did not converge
anova(OTUlm3,OTUlm4)#Mod4 better......................................onlyone that converge

#TSF not significant
UnOtuLm1<- lmer(S.obsU  ~ TimeSinceFire + (1|Plot), data = Unburned, REML=TRUE);summary(UnOtuLm1)#INTcp & 11 sign
UnOTUav<-anova(UnOtuLm1);UnOTUav#p=0.001865

#Post-hoc analysis...................................................
#** Not needed since the reg is not significant 
detach(Unburned)



#********************************************************************************************************************----
#EXPORT RESULTS FROM KRUSKAL WALLIS--------------------------------------------------------------------------------------
#********************************************************************************************************************----

#Treatment .......................................................................................................
capture.output(summary(OTUlm1),file = "Analysis/Diversity/Fungi/Tables/OTU-Trt-lmer.csv")
capture.output(OTUav,file = "Analysis/Diversity/Fungi/Tables/OTU-Trt-anova.csv")

capture.output(summary(OTUlmTrtTSF1),file = "Analysis/Diversity/Fungi/Tables/OTU-Trt-TSF-lmer.csv")
capture.output(OTUavTrtTSF,file = "Analysis/Diversity/Fungi/Tables/OTU-Trt-TSF-anova.csv")


#Burned..............................................................................................................
capture.output(summary(BurnOtuLm1),file = "Analysis/Diversity/Fungi/Tables/Burned-TSF-lmer.csv")
capture.output(BurnOTUav,file = "Analysis/Diversity/Fungi/Tables/Burned-TSF-anova.csv")


#Unburned.............................................................................................................
capture.output(summary(UnOtuLm1),file = "Analysis/Diversity/Fungi/Tables/Unburned-TSF-lmer.csv")
capture.output(UnOTUav,file = "Analysis/Diversity/Fungi/Tables/Unburned-TSF-anova.csv")




