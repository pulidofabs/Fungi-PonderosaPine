#Sept 6, 2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY-----------------------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")

#Load Libraries----------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(data.table)# to calculate percent change


#Load Data------------------------------------------------------------------------------------------------------------
MetaRare<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv", na.strings = "N/A", header = TRUE) #contains species richness metrics
EcMRare3<-read.csv("Analysis/ASV-Tables/Ectomycorrhizal/RareEcM.csv", row.names = 1, check.names = FALSE)

##
##
#*****************************************************************************************************************----
# QUALITY CONTROL )---------------------------------------------------------------------------------------------------
#*****************************************************************************************************************----
# Convert to factor.................................................
MetaRare$FireYear <- as.factor(MetaRare$FireYear)
MetaRare$TimeSinceFire <- as.factor(MetaRare$TimeSinceFire)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRare)
Unburned<-MetaRare[which(Treatment== "Unburned"), ]
head(Unburned);dim(Unburned)#105x33

# * * BURNED PLOTS ------------------------------
Burned<-MetaRare[which(Treatment== "Burned"), ]
head(Burned);dim(Burned)#87x33
detach(MetaRare)

# ALPHA DIVERSITY- MEAN CALCULATIONS--------------------------------------------------------------------------------------------
MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment <- try(relevel(MetaRare$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison


#
#*****************************************************************************************************************************----
#DESRIPTIVE STATISTICS----------------------------------------------------------------------------------------------------------
#*****************************************************************************************************************************----

TrtStats<-MetaRare %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TrtStats # B=8.413793, U=15.514286

TSFstats<-MetaRare %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment, TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstats 


#******************************************************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS------------------------------------------------------------------------
#******************************************************************************************************************************----
EcMTrt <-  MetaRare%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);EcMTrt 

EcMTsfBurned <- Burned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T), 
            sd = sd(TC, na.rm=TRUE)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);EcMTsfBurned

EcMTsfUnburned <- Unburned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);EcMTsfUnburned



# EXPORT FILES----------------------------------------------------------------------------------------------------------------
dir.create("Analysis/Diversity/Ectomycorrhizal/Tables")
dir.create("Analysis/Diversity/Ectomycorrhizal/Tables/DescriptiveStats")

capture.output(TrtStats, file="Analysis/Diversity/Ectomycorrhizal/Tables/DescriptiveStats/Alpha-DesciptiveStats-Trt.csv") 
capture.output(TSFstats, file="Analysis/Diversity/Ectomycorrhizal/Tables/DescriptiveStats/Alpha-DesciptiveStats-TSF.csv") 
capture.output(PerChange, file="Analysis/Diversity/Ectomycorrhizal/Tables/DescriptiveStats/PercentChangeYr-Burned.csv") 




#******************************************************************************************************************************----
# CHECK DATA FOR NORMALITY---------------------------------------------------------------------------------------------------------
#******************************************************************************************************************************----
#Check data normality to decide between anova and KW test

par(mfrow=c(1,3))
hist(MetaRare$S.obs)
hist(Burned$S.obs)
hist(Unburned$S.obs)

#Shapiro test ..................................................
shapiro.test(MetaRare$S.obs)#2.429e-05
shapiro.test(Burned$S.obs)#0.004059
shapiro.test(Unburned$S.obs)#0.389------------#normal ANOVA


#Transform data for normality ...................................
OTUsq<-sqrt(MetaRare$S.obs)
      hist(OTUsq);shapiro.test(OTUsq)#0.02753

BurnOTU<-sqrt(Burned$S.obs)
hist(BurnOTU);shapiro.test(BurnOTU)
      
#******** Data was normalized, so will use lm and anova to test for significance    

#
#*******************************************************************************************************************************----
# PLOTS SPECIES RICHNESS (ASV's)----------------------------------------------------------------------------------------------------
#*******************************************************************************************************************************----

# TREATMENT ----------------------------------------------------------------------------------------------
attach(MetaRare)
SppTrtEcM<-ggplot(MetaRare, aes(x=Treatment, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  ylim(0,30)+
  labs(x = "Sites", y= "Species Richness (ASV's)")+ theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f")); SppTrtEcM

SppTrt1EcM<-SppTrtEcM + stat_compare_means(method = "anova", label.y = 30,); SppTrt1EcM

# SPECIES RICHNESS -- SITES -----------------------------------------------------------------------------
SppSiteEcM<-ggplot(MetaRare, aes(x=Site, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+ 
  ylim(0,30)+
  labs(x = "Sites", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#45877f","#a2673f"));SppSiteEcM

#SppSites<-SppSite+stat_compare_means(method = "kruskal", label.y=30, show.legend = FALSE);SppSites

# SPECIES RICHNESS - FIRE YEAR (looks same as old one)-----------------------------------------------------
SppTrtFYecm<- ggplot(MetaRare, aes(x=FireYear, y=S.obs, group=Treatment, col=Treatment))+
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
  theme(legend.position="bottom");SppTrtFYecm


SppTrtFYecm2<-SppTrtFYecm+stat_compare_means(method = "anova", label.y=20, label =  "p.signif", 
                  show.legend = FALSE, size = 10);SppTrtFYecm2


#SPECIES RICHNESS TIMES SINCE FIRE-----------------------------------------------------------------------------------
SppTrtTSFecm<- ggplot(MetaRare, aes(x=TimeSinceFire, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#45877f","#a2673f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Time Since Fire (Years)", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");SppTrtTSFecm


SppTrtTSFecm2<-SppTrtTSFecm+stat_compare_means(method = "anova", label.y=20, 
                    label =  "p.signif", show.legend = FALSE);SppTrtTSFecm2


# BURNED SPECIES RICHNESS W TIME-SINCE-FIRE --------------------------------------------
TSFecmB<-ggplot(Burned, aes(x=TimeSinceFire, y=S.obs, group=1)) +
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
  theme(legend.position="bottom"); TSFecmB

# UNBURNED SPECIES RICHNESS W TIME-SINCE-FIRE ------------------------------------------
TSFecmUn<-ggplot(Unburned, aes(x=TimeSinceFire, y=S.obs, group=1)) +
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
  theme(legend.position="bottom"); TSFecmUn




#
#
#**************************************************************************************************----
#EXPORT PLOTS -----------------------------------------------------------------------------------------
#**************************************************************************************************----
dir.create("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha")


pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/AsV_Treatment.pdf", height=6, width=8)
SppTrt1EcM
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/AsV_Trt-Site.pdf", height=6, width=8)
SppSiteEcM
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/AsV_Trt-FireYear.pdf", height=6, width=8)
SppTrtFYecm2
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/AsV_Trt-TSF.pdf", height=6, width=8)
SppTrtTSFecm2
dev.off()


#*******************************************************************************************************************----
# -------------- CREATE PANELS FOR PUBLISHING --------------------------------------------------------------------------
#*******************************************************************************************************************----
# Ran codes individually for each guild and then used code below to create the merged panels

linears<-ggarrange(SppTrtFYecm2,SimpTrtTSFsap2, labels = "auto", common.legend = TRUE, legend = "bottom");linears



#Export graphs----------------------------------------------------------------------------------
pdf("Analysis/Diversity/PublicationGraphs/LinearSppRichness.pdf", height=8, width=12)
linears
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

#First check the level of nestedness to see which one to use .............................................
# Looks like OTUlm1 is the model to use 
OTUlm1<- lmer(OTUsq ~ Treatment + (1|Site/Plot/Subplot), data = MetaRare, REML=TRUE);summary(OTUlm1)#does not converge
OTUlm2<- lmer(OTUsq ~ Treatment + (1|Site/Plot), data = MetaRare, REML=TRUE);summary(OTUlm2)#Model does not converge
OTUlm3<- lmer(OTUsq  ~ Treatment + (1|Site), data = MetaRare, REML=TRUE);summary(OTUlm3)#...........best model
OTUlm4<- lmer(OTUsq  ~ Treatment + (1|Plot), data = MetaRare, REML=TRUE);summary(OTUlm4)

anova(OTUlm2,OTUlm3)#Mod3 better
anova(OTUlm3,OTUlm4)#Mod3 better

#Model gave the same results so will proceed with Site level nestedness for all analysis


##
####
#*************************************************************************************************************----
#-------------------  TEST FOR SIGNIFICANCE OF THE DATA ----------------------------------------------------------
#*************************************************************************************************************----
#Effect of treatment on species richness (using transform species richness for significance).........................
OTUlmTrtTSF1<- lmer(OTUsq ~ Treatment * TimeSinceFire + (1|Site), data = MetaRare, REML=TRUE);summary(OTUlmTrtTSF1)
OTUavTrtTSF<-anova(OTUlmTrtTSF1);OTUavTrtTSF#Gives us the same pvalue 
#0.0065 ** 

#Final nestedness test .................................................. ..........................
OTUlm1<- lmer(OTUsq ~ Treatment + (1|Site), data = MetaRare, REML=TRUE);summary(OTUlm1)
OTUav<-anova(OTUlm1);OTUav#Gives us the same pvalue 

#Test TimeSinceFire
OTUlm1TSF<- lmer(OTUsq ~ TimeSinceFire + (1|Site), data = MetaRare, REML=TRUE);summary(OTUlm1TSF)
OTUavTSF<-anova(OTUlm1TSF);OTUavTSF

#Post-Hoc Analysis............................................................


#TIME-SINCE-FIRE BURNED DATA -------------------------------------------------------------------------------------
attach(Burned)
OTUlmB1<- lmer(BurnOTU ~ TimeSinceFire + (1|Site), data = Burned, REML=TRUE);summary(OTUlmB1)
OTUlmB2<- lmer(BurnOTU ~ TimeSinceFire + (1|Plot), data = Burned, REML=TRUE);summary(OTUlmB2)
OTUlmB3<- lmer(BurnOTU ~ TimeSinceFire + (1|Site/Plot), data = Burned, REML=TRUE);summary(OTUlmB3)#no conver
AIC(OTUlmB1,OTUlmB2)#otuLM1 better 

BurnOtuLm1<- lmer(BurnOTU ~ TimeSinceFire + (1|Site), data = Burned, REML=TRUE);summary(BurnOtuLm1)#pr>1
BurnOTUav<-anova(BurnOtuLm1);BurnOTUav#Pvalue = 1, not significant 

#Post-hoc analysis...................................................
#** Not needed since the reg is not significant 
detach(Burned)



#TIM-SINCE-FIRE UNBURNED DATA ------------------------------------------------------------------------------------
attach(Unburned)
OTUlmUn1<- lmer(S.obs ~ TimeSinceFire + (1|Site), data = Unburned, REML=TRUE);summary(OTUlmUn1)
OTUlmUn2<- lmer(S.obs ~ TimeSinceFire + (1|Plot), data = Unburned, REML=TRUE);summary(OTUlmUn2)
OTUlmUn3<- lmer(S.obs ~ TimeSinceFire + (1|Site/Plot), data = Unburned, REML=TRUE);summary(OTUlmUn3)#no conver
AIC(OTUlmUn1,OTUlmUn2)#Same AIC

UnOtuLm1<- lmer(S.obs ~ TimeSinceFire + (1|Site), data = Unburned, REML=TRUE);summary(UnOtuLm1)
UnOTUav<-anova(UnOtuLm1);UnOTUav#p=0.3349, Not significant 

#Post-hoc analysis...................................................
#** Not needed since the reg is not significant 
detach(Unburned)



#********************************************************************************************************************----
#EXPORT RESULTS FROM KRUSKAL WALLIS--------------------------------------------------------------------------------------
#********************************************************************************************************************----

#Treatment .......................................................................................................
capture.output(summary(OTUlm1),file = "Analysis/Diversity/Ectomycorrhizal/Tables/OTU-Trt-lmer.csv")
capture.output(OTUav,file = "Analysis/Diversity/Ectomycorrhizal/Tables/OTU-Trt-anova.csv")

capture.output(summary(OTUlmTrtTSF1),file = "Analysis/Diversity/Ectomycorrhizal/Tables/OTU-Trt-TSF-lmer.csv")
capture.output(OTUavTrtTSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/OTU-Trt-TSF-anova.csv")


#Burned..............................................................................................................
capture.output(summary(BurnOtuLm1),file = "Analysis/Diversity/Ectomycorrhizal/Tables/Burned-TSF-lmer.csv")
capture.output(BurnOTUav,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Burned-TSF-anova.csv")


#Unburned.............................................................................................................
capture.output(summary(UnOtuLm1),file = "Analysis/Diversity/Ectomycorrhizal/Tables/Unburned-TSF-lmer.csv")
capture.output(UnOTUav,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Unburned-TSF-anova.csv")



#
#
#*****************************************************************************************************************----
#--------------------    SAME CALCULATIONS AS ABOVE USING SIMPSONS DIVERSITY INDEX -----------------------------------
#*****************************************************************************************************************----
attach(MetaRare)

par(mfrow=c(1,3))
hist(MetaRare$simpson);shapiro.test(MetaRare$simpson)
hist(Burned$simpson); shapiro.test(Burned$simpson)
hist(Unburned$simpson); shapiro.test(Unburned$simpson)

#Transformations
Simplog<-log(MetaRare$simpson); shapiro.test(Simplog)
SimpBlog<-log(Burned$simpson);shapiro.test(SimpBlog)
SimpUNlog<-log(Unburned$simpson); shapiro.test(SimpUNlog)



# TREATMENT -------------------------------------------------------------------------------------------
SppTrtSimpEcM<-ggplot(MetaRare, aes(x=Treatment, y=simpson)) +   
              geom_boxplot(aes(fill=Treatment))+ 
              theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,15)+  
              labs(x = "Sites", y= "Simpson's Index")+ theme(legend.position="bottom")+
              scale_fill_manual(values=c("#45877f","#a2673f")); SppTrtSimpEcM

SppTrtSimpEcM1<-SppTrtSimpEcM + stat_compare_means(method = "kruskal", label.y = 15,); SppTrtSimpEcM1

# SPECIES RICHNESS -- SITES -----------------------------------------------------------------------------
SppSiteSimpEcM<-ggplot(MetaRare, aes(x=Site, y=simpson)) +   
                geom_boxplot(aes(fill=Treatment))+ 
                theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,15)+   
                labs(x = "Sites", y= "Simpson's Index")+theme(legend.position="bottom")+
                scale_fill_manual(values=c("#45877f","#a2673f"));SppSiteSimpEcM


# SPECIES RICHNESS - FIRE YEAR (looks same as old one)-----------------------------------------------------------------------------------------
SppTrtFYSimpEcM<- ggplot(MetaRare, aes(x=FireYear, y=simpson, group=Treatment, col=Treatment))+
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
                  theme(legend.position="bottom");SppTrtFYSimpEcM

SppTrtFYsimpEcM2<-SppTrtFYSimpEcM+stat_compare_means(method = "kruskal", label.y=9, label =  "p.signif", 
                           show.legend = FALSE);SppTrtFYsimpEcM2


#SPECIES RICHNESS TIMES SINCE FIRE-----------------------------------------------------------------------------------
SppTrtTSFsimpEcM<- ggplot(MetaRare, aes(x=TimeSinceFire, y=simpson, group=Treatment, col=Treatment))+
                  stat_summary(fun=mean,geom="line", size=1.3)+
                  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
                               alpha=0.7,position = position_dodge(0.01))+
                  scale_color_manual(values=c("#45877f","#a2673f")) + 
                  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=.7),
                     text = element_text(size=23),
                     axis.text.y = element_text(size=19,hjust=0.5),
                     axis.text.x = element_text(size=19))+
                  labs(x = "Time Since Fire (Years)", y= "Simpson Index")+
                  theme(legend.position="bottom");SppTrtTSFsimpEcM

SppTrtTSFsimpEcM2<-SppTrtTSFsimpEcM+stat_compare_means(method = "kruskal", label.y=10, 
                   label =  "p.signif", show.legend = FALSE);SppTrtTSFsimpEcM2

#
#
#**************************************************************************************************************----
#EXPORT PLOTS -----------------------------------------------------------------------------------------------
#**************************************************************************************************************----
dir.create("Analysis/Diversity/Ectomycorrhizal/Graphs")

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/Simpson_Alpha_Treatment.pdf", height=6, width=8)
SppTrtSimpEcM1
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/Simpson_Alpha_Trt-Site.pdf", height=6, width=8)
SppSiteSimpEcM
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/Simpson_Alpha_Trt-FireYear.pdf", height=6, width=8)
SppTrtFYsimpEcM2
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/Alpha/Simpson_Alpha_Trt-TSF.pdf", height=6, width=8)
SppTrtTSFsimpEcM2
dev.off()


#
#
#*************************************************************************************************************----
#TEST SIGNIFICANCE -----------------------------------------------------------------------------------------------
#*************************************************************************************************************----
# TREATMENT -------------Kruskal Wallis-treatment-all data (no subsetting)
attach(MetaRare)
SimpkwTrt<-kruskal.test(simpson~Treatment, MetaRare); SimpkwTrt #p<4.31e-07
detach(MetaRare)

#TIM-SINCE-FIRE UNBURNED DATA ------------------------------------------------------------------------------------
attach(Unburned)
Un.Simp.TSF<-kruskal.test(simpson~TimeSinceFire, Unburned);Un.Simp.TSF #p=1.397e-05
Pair.Simp.UN.TSF<-pairwise.wilcox.test(simpson, TimeSinceFire, p.adjust.method = "bonf", data=Unburned)
Pair.Simp.UN.TSF #Significance 3-5; 3-11
detach(Unburned)

#TIME-SINCE-FIRE BURNED DATA -------------------------------------------------------------------------------------
attach(Burned)
B.Simp.TSF<-kruskal.test(simpson~TimeSinceFire, Burned);B.Simp.TSF #p=0.03412 Significant
Pair.Simp.B.TSF<-pairwise.wilcox.test(simpson, TimeSinceFire, p.adjust.method = "bonf", data=Burned);
Pair.Simp.B.TSF # 2-11 only
detach(Burned)


#**************************************************************************************************************----
#EXPORT RESULTS FROM KRUSKAL WALLIS--------------------------------------------------------------------------------
#**************************************************************************************************************----
capture.output(SimpkwTrt,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Simpson-Treatment-KW.CSV")

#Unburned
capture.output(Un.Simp.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Simpson-TSF-KW-Unburned.CSV")
capture.output(Pair.Simp.UN.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Simpson-TSF-pairwise-Unburned.CSV")

#Burned
capture.output(B.Simp.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Simpson-TSF-KW-Burned.CSV")
capture.output(Pair.Simp.B.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Sipson-TSF-pairwise-Burned.CSV")













































