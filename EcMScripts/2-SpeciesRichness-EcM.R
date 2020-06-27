#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY-----------------------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

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
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);EcMTrt 

EcMTsfBurned <- Burned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
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

#shapiro test 
shapiro.test(MetaRare$S.obs)#2.429e-05
shapiro.test(Burned$S.obs)#0.004059
shapiro.test(Unburned$S.obs)#0.389------------#normal ANOVA


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
  scale_fill_manual(values=c("#a2673f", "#45877f")); SppTrtEcM

SppTrt1EcM<-SppTrtEcM + stat_compare_means(method = "kruskal", label.y = 30,); SppTrt1EcM

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
  scale_fill_manual(values=c("#a2673f", "#45877f"));SppSiteEcM

#SppSites<-SppSite+stat_compare_means(method = "kruskal", label.y=30, show.legend = FALSE);SppSites

# SPECIES RICHNESS - FIRE YEAR (looks same as old one)-----------------------------------------------------
SppTrtFYecm<- ggplot(MetaRare, aes(x=FireYear, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
               alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#a2673f", "#45877f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Time Since Fire (Years)", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");SppTrtFYecm


SppTrtFYecm2<-SppTrtFYecm+stat_compare_means(method = "kruskal", label.y=20, label =  "p.signif", 
                                       show.legend = FALSE, size = 10);SppTrtFYecm2


#SPECIES RICHNESS TIMES SINCE FIRE-----------------------------------------------------------------------------------
SppTrtTSFecm<- ggplot(MetaRare, aes(x=TimeSinceFire, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1.3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(values=c("#a2673f", "#45877f")) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7),
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
  labs(x = "Time Since Fire (Years)", y= "Species Richness (ASV's)")+
  theme(legend.position="bottom");SppTrtTSFecm


SppTrtTSFecm2<-SppTrtTSFecm+stat_compare_means(method = "kruskal", label.y=20, 
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
attach(MetaRare)
OTUkwTrt<-kruskal.test(S.obs~Treatment, MetaRare); OTUkwTrt #p<2.2e-16
detach(MetaRare)

#TIM-SINCE-FIRE UNBURNED DATA ------------------------------------------------------------------------------------
attach(Unburned)
Un.OTU.TSF<-kruskal.test(S.obs~TimeSinceFire, Unburned);Un.OTU.TSF #p=4.288e-05
Pair.OTU.UN.TSF<-pairwise.wilcox.test(S.obs, TimeSinceFire, p.adjust.method = "bonf", data=Unburned)
Pair.OTU.UN.TSF #Significance 2-5; 3-5
detach(Unburned)

#TIME-SINCE-FIRE BURNED DATA -------------------------------------------------------------------------------------
attach(Burned)
B.OTU.TSF<-kruskal.test(S.obs~TimeSinceFire, Burned);B.OTU.TSF #p=0.09409 NOT SIGNIFICANT
Pair.OTU.B.TSF<-pairwise.wilcox.test(S.obs, TimeSinceFire, p.adjust.method = "bonf", data=Burned);Pair.OTU.B.TSF
detach(Burned)


#**************************************************************************************************************----
#EXPORT RESULTS FROM KRUSKAL WALLIS--------------------------------------------------------------------------------
#**************************************************************************************************************----
capture.output(OTUkwTrt,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Otu-Treatment-KW.CSV")

#Unburned
capture.output(Un.OTU.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Otu-TSF-KW-Unburned.CSV")
capture.output(Pair.OTU.UN.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Otu-TSF-pairwise-Unburned.CSV")

#Burned
capture.output(B.OTU.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Otu-TSF-KW-Burned.CSV")
capture.output(Pair.OTU.B.TSF,file = "Analysis/Diversity/Ectomycorrhizal/Tables/Otu-TSF-pairwise-Burned.CSV")



#
#
#*****************************************************************************************************************----
#--------------------    SAME CALCULATIONS AS ABOVE USING SIMPSONS DIVERSITY INDEX -----------------------------------
#*****************************************************************************************************************----

# TREATMENT -------------------------------------------------------------------------------------------
attach(MetaRare)
SppTrtSimpEcM<-ggplot(MetaRare, aes(x=Treatment, y=simpson)) +   
              geom_boxplot(aes(fill=Treatment))+ 
              theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                   text = element_text(size=23),
                   axis.text.y = element_text(size=19,hjust=0.5),
                   axis.text.x = element_text(size=19))+ ylim(0,15)+  
              labs(x = "Sites", y= "Simpson's Index")+ theme(legend.position="bottom")+
              scale_fill_manual(values=c("#a2673f", "#45877f")); SppTrtSimpEcM

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
                scale_fill_manual(values=c("#a2673f", "#45877f"));SppSiteSimpEcM


# SPECIES RICHNESS - FIRE YEAR (looks same as old one)-----------------------------------------------------------------------------------------
SppTrtFYSimpEcM<- ggplot(MetaRare, aes(x=FireYear, y=simpson, group=Treatment, col=Treatment))+
                  stat_summary(fun=mean,geom="line", size=1.3)+
                  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85,
                               alpha=0.7,position = position_dodge(0.01))+
                  scale_color_manual(values=c("#a2673f", "#45877f")) + 
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
                  scale_color_manual(values=c("#a2673f", "#45877f")) + 
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













































