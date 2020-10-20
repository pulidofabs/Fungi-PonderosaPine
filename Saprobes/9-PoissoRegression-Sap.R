#10-9-2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
#setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")


#Load Libraries----------------------
library(ggpubr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lsmeans)
library(lme4)
library(ochRe)
library(rcompanion)


#Load Data--------------------------------------------------------------------------
MetaRareSap<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.csv",
                    na.strings = "N/A", header = TRUE) 


#*******************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#*******************************************************************************************************----
# Convert to factor------------------------------------------------
MetaRareSap$FireYear <- as.factor(MetaRareSap$FireYear)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRareSap)

# * * UNBURNED PLOTS ------------------------------
UnburnedSap<-MetaRareSap[which(Treatment== "Unburned"), ]
head(UnburnedSap);dim(UnburnedSap)#102x34

# * * BURNED PLOTS ------------------------------
BurnedSap<-MetaRareSap[which(Treatment== "Burned"), ]
head(BurnedSap);dim(BurnedSap)#106x34
detach(MetaRareSap)


#***************************************************************************************************************----
# POISSON REGRESSION -----------------------------------------------------------------------------------------------
#***************************************************************************************************************----
#we assume that the dependent variable is over-dispersed and does not have an excessive number of zeros. 
attach(BurnedSap)
attach(MetaRareSap)

#CHECK ASSUMPTIONS ................................................................
#Check the means and SD of the data to verify that neg binom is right

with(MetaRareSap, tapply(S.obs,TimeSinceFire,  function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))



#Run the null model to veridy level of nestedness ................................................................................
ASVpoissonSap1 <- glmer(S.obs ~ TimeSinceFire + (1|Site/Plot/Subplot), family="poisson",data = BurnedSap);summary(ASVpoissonSap1)# no conver
ASVpoissonSap2 <- glmer(S.obs ~ TimeSinceFire + (1|Site/Plot), family="poisson",data = BurnedSap);summary(ASVpoissonSap2)#no converg
ASVpoissonSap3 <- glmer(S.obs ~ TimeSinceFire + (1|Site), family="poisson",data = BurnedSap);summary(ASVpoissonSap3)
ASVpoissonSap4 <- glmer(S.obs ~ TimeSinceFire + (1|Plot), family="poisson",data = BurnedSap);summary(ASVpoissonSap4)

#Sap1 and Sap2 did not convergence, so will not use, 
anova(ASVpoissonSap3,ASVpoissonSap4)#Sap4 better but no convergence
anova(ASVpoissonSap1,ASVpoissonSap2, ASVpoissonSap3,ASVpoissonSap4)



#RUN THE POISSON MODEL -----------------------------------------------------------------
#Mean is higher than the variance (SD)-- so do not need to perform neg binom, 
#instead we will run a poisson regression 
library(emmeans)

#overall model P-value
ASVpoissonSap <- glm(S.obs ~ TimeSinceFire + (1/Plot), family=poisson, data = BurnedSap);summary(ASVpoissonSap)

#Make into a factor so that we can look at contrast
ASVpoissonSap <- glm(S.obs ~ as.factor(TimeSinceFire) + (1/Plot), family=poisson, data = BurnedSap);summary(ASVpoissonSap)

ASVemm <- emmeans(ASVpoissonSap, "TimeSinceFire"); ASVemm
pairs(ASVemm)
pwpp(ASVemm)






#EXPORT POISSON RESULTS --------------------------------------------------------------------------------------
dir.create("Analysis/MixedModels/Saprobes/Poisson")
capture.output(summary(ASVpoissonSap), file="Analysis/MixedModels/Saprobes/Poisson/RegressionPoisson-Saprobes-glm.csv")


#*****************************************************************************************************************----
#------------------   MODEL VALIDATION AND FINAL PLOT ----------------------------------------------------------------
#*****************************************************************************************************************----
par(mfrow=c(2,2))

#Model validations
plot(ASVpoissonSap)




#PLOT FINAL GRAPH--------------------------------------------------------------------------------------------
MetaRareSap$FireYear<-as.factor(MetaRareSap$TimeSinceFire)

TSFasvSap<-ggplot(MetaRareSap, aes(x=TimeSinceFire, y=S.obs))  +
  geom_smooth(method = "glm", se = TRUE, method.args = list(family = "poisson"),
              aes(color=Treatment))+
  geom_point(aes(fill=factor(FireYear)), size=4, shape=21, stroke=0) +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
           sep = "~~~~")), formula = y~x, label.x = 1, label.y = 200, size=6)+
  scale_fill_ochre(palette="olsen_seq", discrete = TRUE)+
  scale_colour_manual(values=c("#a2673f","#45877f")) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        axis.text.y = element_text(size =25, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 25,  hjust=1, colour = "black",
                 margin = margin(t = 0, r = 0, b = 20, l = 0)), 
        axis.title.y = element_text(size = 25, colour = "black",
                 margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.title = element_text(size=28),
        legend.text=element_text(size = 22),
        legend.position = "bottom")+ 
  labs(fill="Time Since Fire (Years)", colour="Treatment")+
  labs(x = "Time Since Fire (Years)",y="Species Richness");TSFasvSap

  
TSFasvSap1<-TSFasvSap + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  guides(color = guide_legend(override.aes = list(size=4))) ;TSFasvSap1




#
#
#*********************************************************************************************************----
# --- EXPORT GRAPHS ------------------------------------------------------------------------------------------
#*********************************************************************************************************----

pdf("Analysis/MixedModels/Saprobes/Poisson/ASV-Sap-Poisson-TSF-Trt.pdf", height=6, width=7)
TSFasvSap1
dev.off()



#DO NOT CLOSE OR CLEAR ENVIRONMENT NEED RESULTS TO CREATE PANELS FOR PUBLICATION










