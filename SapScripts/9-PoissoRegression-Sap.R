#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

#Load Libraries----------------------
library(ggpubr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lsmeans)


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

#CHECK ASSUMPTIONS ................................................................
#Check the means and SD of the data to verify that neg binom is right

with(BurnedSap, tapply(S.obs,TimeSinceFire,  function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))



#RUN THE POISSON MODEL -----------------------------------------------------------------
#Mean is higher than the variance (SD)-- so do not need to perform neg binom, 
#instead we will run a poisson regression 
ASVpoissonSap <- glm(S.obs ~ TimeSinceFire + (1/Plot), family="poisson",data = BurnedSap);ASVpoissonSap
summary(ASVpoissonSap) #0.000266
AovSap<-anova(ASVpoissonSap); summary(AovSap)#TSF is sigificant
nagelkerke(ASVpoissonSap)

#EXPORT POISSON RESULTS --------------------------------------------------------------------------------------
dir.create("Analysis/MixedModels/Saprobes/Poisson")
capture.output(summary(ASVpoissonSap), file="Analysis/MixedModels/Saprobes/Poisson/RegressionPoisson-Saprobes.csv")


#*****************************************************************************************************************----
#------------------   MODEL VALIDATION AND FINAL PLOT ----------------------------------------------------------------
#*****************************************************************************************************************----
par(mfrow=c(2,2))
plot(ASVpoissonSap)


#Modelvalidation 
par(mfrow=c(1,1))
F1Sap <- fitted(ASVpoissonSap) #command fitted gives already e^model
E1Sap <- resid(ASVpoissonSap, type = "pearson")
plot(x = F1Sap, 
     y = E1Sap,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")




#PLOT FINAL GRAPH--------------------------------------------------------------------------------------------
BurnedSap$FireYear<-as.factor(FireYear)

TSFasvSap<-ggplot(BurnedSap, aes(TimeSinceFire, log(S.obs)))+
              geom_point(aes(col= FireYear, size=.7))+ 
              scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
              geom_smooth(method = "glm", se = TRUE,  
                  method.args = list(family = "poisson"))+
               stat_cor(label.x = 9, label.y = 5.5) +
              stat_regline_equation(label.x = 9, label.y = 5.4)+
              theme_bw()+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=.7), 
                    text = element_text(size=23),
                    axis.text.y = element_text(size=19,hjust=0.5),
                    axis.text.x = element_text(size=19))+
              labs(x = "Time Since Fire (years)", y= "")+ 
              theme(legend.position="bottom");TSFasvSap


TSFasvSap1<-TSFasvSap + scale_x_discrete(breaks = scales::pretty_breaks(n = 10))+
  guides(color = guide_legend(override.aes = list(size=4)));TSFasvSap1




#
#
#*********************************************************************************************************----
# --- EXPORT GRAPHS ------------------------------------------------------------------------------------------
#*********************************************************************************************************----

pdf("Analysis/MixedModels/Saprobes/Poisson/ASV-Sap-Poisson-TSF.pdf", height=6, width=7)
TSFasvSap1
dev.off()



#DO NOT CLOSE OR CLEAR ENVIRONMENT NEED RESULTS TO CREATE PANELS FOR PUBLICATION










