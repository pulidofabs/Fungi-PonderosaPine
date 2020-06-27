#jUNE 24, 2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

#Load libraries---------
library(ggplot2)
library(ggpubr)
library(dplyr)



#Load Data--------------------------------------------------------------------------
MetaRare<- read.csv("Analysis/Metadata/Fungi/MetaRareCore.csv", na.strings = "N/A", header = TRUE) #Rarefied metadata


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

detach(MetaRare)

#
#
#***********************************************************************************************----
# ALPHA DIVERSITY- MEAN CALCULATIONS----------------------------------------------------------------
#***********************************************************************************************----

MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment <- try(relevel(MetaRare$Treatment , "Unburned"))



TrtStats<-MetaRare %>%
  filter(!is.na(ObsOtu)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),Avg_ASV = mean(ObsOtu, na.rm = TRUE), 
            Variance = var(ObsOtu, na.rm = TRUE), sd = sd(ObsOtu), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TrtStats # B=203.6827, U=232.7925

TSFstats<-MetaRare %>%
  filter(!is.na(ObsOtu)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(ObsOtu, na.rm = TRUE), 
            Variance = var(ObsOtu, na.rm = TRUE), sd = sd(ObsOtu), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstats # 2-205-2963, 3=204.7925, 5= 216.6604, 11=248.7200 


TSFstatsB<-Burned %>%
  filter(!is.na(ObsOtu)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(ObsOtu, na.rm = TRUE), 
            Variance = var(ObsOtu, na.rm = TRUE), sd = sd(ObsOtu), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 1) %>% as.data.frame();TSFstatsB # 2-202.5185, 3=228.6538, 5= 199.6923, 11=183.1200


TSFstatsUn<-Unburned %>%
  filter(!is.na(ObsOtu)) %>%
  group_by(TimeSinceFire) %>%
  summarise(n_obs = n(),Avg_ASV = mean(ObsOtu, na.rm = TRUE), 
            Variance = var(ObsOtu, na.rm = TRUE), sd = sd(ObsOtu), se =sd / sqrt(n_obs)) %>%
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
  summarise(mean = mean(ObsOtu,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);FungiTrt #-12.5

FungiTsfBurned <- Burned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(ObsOtu,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);FungiTsfBurned

FungiTsfUnburned <- Unburned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(ObsOtu,na.rm = T)) %>% 
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
par(mfrow=c(1,3))
hist(MetaRare$ObsOtu)
hist(Burned$ObsOtu)
hist(Unburned$ObsOtu)

#shapiro test 
shapiro.test(MetaRare$ObsOtu)#0.0002785
shapiro.test(Burned$ObsOtu)#0.0163
shapiro.test(Unburned$ObsOtu)#0.04436

#tranform to see if it makes a difference, if yes = ANOVA if no WILCOXON TEST
ObsOTU<-log(MetaRare$ObsOtu)
ObsOTUB<-log(Burned$ObsOtu)
ObsOTUU<-log(Unburned$ObsOtu)

shapiro.test(ObsOTU)#0.578
shapiro.test(ObsOTUB)#0.6702
shapiro.test(ObsOTUU)#0.2061


#
#
#-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*_*_*_*
# PLOTS SPECIES RICHNESS (ESV's)---------------------------------------------------------------------
#-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*_*_*_*

# TREATMENT -------------------------------------------------------------------------------------------
SppTrt<-ggboxplot(MetaRare, x="Treatment", y="ObsOTU") +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw() +  ylab("log (Species Richness (ASV's))") +  
  theme(legend.position = "bottom",text = element_text(size=20),
        axis.title.x=element_blank(), axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18))+
  scale_fill_manual(values=c( "#45877f","#a2673f")); SppTrt

SppTrt1<-SppTrt + stat_compare_means(method = "anova", label.y = 6.2); SppTrt1

# SPECIES RICHNESS -- SITES --------------------------------------------------------------------------------------
SppTSFBox<-ggplot(MetaRare, aes(x=TimeSinceFire, y=ObsOTU)) +   
           geom_boxplot(aes(fill=Treatment))+ 
            theme_bw()+
            theme(legend.position = "bottom", 
                  text = element_text(size=20),
                  axis.text.y = element_text(size=20), 
                  axis.text.x = element_text(size=20))+
            scale_fill_manual(label=c("Unburned","Burned"), values=c("#45877f","#a2673f")) +
            labs(x = "Time Since Fire (Years)", y= "Log (Species Richness (ASV's))");SppTSFBox

#SppFireY<-SppFY + stat_compare_means(method = "anova", paired=TRUE,label.y = 6.05);SppFireY



# SPECIES RICHNESS - TSF -----------------------------------------------------------------------------------------
SppTSF<-ggplot(MetaRare, aes(x=TimeSinceFire, y=ObsOTU, group=Treatment, col=Treatment))+
              stat_summary(fun=mean,geom="line", size=.85)+
              stat_summary(fun.data = mean_se,geom = "errorbar", size=.5, alpha=0.7,position = position_dodge(0.01))+
              theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                    text = element_text(size=24),
                    axis.text.y = element_text(size=22,hjust=0.5),
                    axis.text.x = element_text(size=22))+ scale_color_manual(values=c("#45877f","#a2673f")) +   
                    labs(x = "Time-Since-Fire (days)", y= "Species Richness (ASV's)")+
                theme(legend.position="bottom");SppTSF
SppTSF1<-SppTSF + stat_compare_means(method = "anova", label.y = 5.8, label = "p.signif", size=6,
                                     show.legend=FALSE);SppTSF1



# EXPORT PLOTS ---------------------------------------------------------------------------------------------------------


pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-Treatment.pdf", height=6, width=8)
SppTrt1
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-TSF-Boxplot.pdf", height=6, width=8)
SppTSFBox
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/SpeciesRichness-TSF-Line.pdf", height=8, width=10)
SppTSF1
dev.off()


#---------CREATE A PANEL
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



#**************************************************************************************************
# - TEST FOR SIGNIFICANCE 
#**************************************************************************************************

par(mfrow=c(1,3))
hist(MetaRare$ObsOtu)
hist(MetaRare$simpson)

shapiro.test(log(Burned$ObsOtu))# log works
shapiro.test(log(MetaRare$ObsOtu))# log makes it normal--anova
shapiro.test((MetaRare$simpson))# cannot normalize- KW


OTU<-log(MetaRare$ObsOtu)
OtuB<-log(Burned$ObsOtu)

#Species richness -----------------------------------------------
OtuTrtlm <- lm(OTU ~ Treatment, data = MetaRare);OtuTrtlm
OtuAV <- aov(OtuTrtlm);summary(OtuAV) #0.00133 * * 


OtuTsFlm <- lm(OtuB ~ TimeSinceFire, data = Burned);OtuTsFlm
OtuAVtsf <- aov(OtuTsFlm);summary(OtuAVtsf) #0.0227 * 
tukeyOtutsf <- TukeyHSD(OtuAVtsf); tukeyOtutsf

#Export results-------------------------------------------------

capture.output(summary(OtuAV), file="Analysis/Diversity/Fungi/Tables/Fungal-SpeciesRichness-Ttr-Anova.csv")
capture.output(summary(OtuAVtsf), file="Analysis/Diversity/Fungi/Tables/Fungal-SpeciesRichness-TSF-Anova.csv")
capture.output(tukeyOtutsf, file="Analysis/Diversity/Fungi/Tables/Fungal-SpeciesRichness-Tukey.csv")





# TREATMENT -------------Kruskal Wallis-treatment-all data (no subsetting)
attach(MetaRare)
SimpsonTrt<-kruskal.test(simpson~Treatment, MetaRare);SimpsonTrt#7.055e-09


#TIM-SINCE-FIRE UNBURNED DATA ------------------------------------------------------------------------------------
attach(Burned)
simpsonTSF<-kruskal.test(simpson~TimeSinceFire, Burned);simpsonTSF #Not sign
SimpsonTSFpair<-pairwise.wilcox.test(simpson, TimeSinceFire, p.adjust.method = "bonf", data=Burned); SimpsonTSFpair
detach(Burned)



#****************************************************************************************---
#- Percent change 
#***************************************************************************************---

OtuTrt <-  MetaRare%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(ObsOtu,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);OtuTrt #12.4 decrease 


Unburned%>% 
  group_by(TimeSinceFire) %>%
summarise(mean = mean(ObsOtu,na.rm = T))

OtuTsf <- Burned%>% 
  group_by(TimeSinceFire) %>% 
  summarise(mean = mean(ObsOtu,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100);OtuTsf 








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












