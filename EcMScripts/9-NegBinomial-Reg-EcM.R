#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

#LOAD LIBRARIES .................................................................
library(ggplot2)
library(tidyr)
library(dplyr)
library(ochRe)
library(MASS)
library(rcompanion)# for pseudo R


#LOAD DATA....................................................................... 
MetaRareEcM<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv",
                    na.strings = "N/A", header = TRUE) 


#***************************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#***************************************************************************************************************----
# Convert to factor------------------------------------------------
MetaRareEcM$FireYear<- as.factor(MetaRareEcM$FireYear)



# SUBSET DATA -----------------------------------------------------
attach(MetaRareEcM)

# * * UNBURNED PLOTS -----------------------------------------
UnburnedEcM<-MetaRareEcM[which(Treatment== "Unburned"), ]
head(UnburnedEcM);dim(UnburnedEcM)#105x33

# * * BURNED PLOTS -------------------------------------------
BurnedEcM<-MetaRareEcM[which(Treatment== "Burned"), ]
head(BurnedEcM);dim(BurnedEcM)#87x33
detach(MetaRareEcM)


#*******************************************************************************************************----
# NEGATIVE BINOMIAL REGRESSION -----------------------------------------------------------------------------
#*******************************************************************************************************----
#we assume that the dependent variable is over-dispersed and does not have an excessive number of zeros. 
attach(BurnedEcM)

levels(BurnedEcM$TimeSinceFire)


#Check the means and variance of the data to verify that neg binom is right
#variance is higher that mean= negative binomial 
df2 <- BurnedEcM %>%
  group_by(TimeSinceFire) %>%
  summarise(mcount = mean(S.obs),
            varcount = var(S.obs)) 


#RUN REGRESSION -------------------------------------------------------------------
library(MuMIn) #to get pseudo R
library(lme4)


ASVneg <- glmer.nb(S.obs ~ TimeSinceFire + (1|Plot),data = BurnedEcM);ASVneg
Res<-summary(ASVneg);Res #0.00628 y = 2.29663X - 0.03492


ASVnull <- glmer.nb(S.obs ~ 1 +( 1 | Plot), data = BurnedEcM)
anova(ASVneg, ASVnull)
AovRes<-anova(ASVneg); AovRes

#CALCULATE R2 FOR THE MODEL -------------------------------------------------------
#R2M = marginal: for fixed effect  and  R2c= conditional for entire model 

R<-r.squaredGLMM(ASVneg, ASVnull);R


#PREDICTED GRAPHS------------------------------------------------------------------ 
library(ggeffects)
pr1 <- ggpredict(ASVneg, c("TimeSinceFire"))
plot(pr1)


#EXPORT RESULTS -----------------------------------------------------------------------------------------
capture.output(Res, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto.csv")
capture.output(AovRes, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto-ANOVA.csv")
capture.output(R, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto-Rsquare.csv")
capture.output(summary(ASVneg), file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/ASV-NegBinom-ModelResults.csv")

##
##
#*********************************************************************************************************----
#--------------- MODEL VALIDATION ----------------------------------------------------------------------------
#*********************************************************************************************************----
BurnedEcM$FireYear<-as.factor(BurnedEcM$FireYear)

#Get the confidence interval levels
(est <- cbind(Estimate = coef(ASVneg), confint(ASVneg)))
exp(est)#Incidence ratio



#Diagnostic plots ..................................................................
par(mfrow=c(1,1))
plot(ASVneg)


#Modelvalidation ...................................................................
par(mfrow=c(1,1))
FittedEcM <- fitted(ASVneg) #command fitted gives already e^model
ResEcM <- resid(ASVneg, type = "pearson")
plot(x = FittedEcM, 
     y = ResEcM,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")



#*******************************************************************************************----
#CREATE GRAPHS ---------------------------------------------------------------------------------
#*******************************************************************************************----
#install.packages("ggpmisc")
library(ggpmisc)
library(MASS)

TSFasvEcM<-ggplot(BurnedEcM, aes(TimeSinceFire, log(S.obs)))+
  geom_point(aes(col= FireYear, size=.5))+ 
  scale_colour_ochre(palette="olsen_seq", discrete = TRUE)+
  geom_smooth(method = MASS::glm.nb, se = TRUE, size=1.5)+
  stat_cor(label.x = 9, label.y = 3) +
  stat_regline_equation(label.x = 9, label.y = 2.85)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=23),
        axis.text.y = element_text(size=19,hjust=0.5),
        axis.text.x = element_text(size=19))+
        labs(x = "Time Since Fire (years)", 
             y= "log(Species Richness (ASV's))")+
  theme(legend.position="none");TSFasvEcM

TSFasvEcM1<-TSFasvEcM + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  guides(color = guide_legend(override.aes = list(size=4))) ;TSFasvEcM1




#
#*********************************************************************************************************----
# --- EXPORT GRAPHS ------------------------------------------------------------------------------------------
#*********************************************************************************************************----


pdf("Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/ASV-NB-TSF-EcM-new.pdf", height=6, width=7)
TSFasvEcM1
dev.off()





#***********************************************************************************************----
# CREATE PANEL OF REGRESSION FOR ECM AND SAPROBES --------------------------------------------------
# **********************************************************************************************----

#Individual scripts for regresssion were first ran for EcM and Saprobe individually, 
#data was kept in environment and then this code was ran to create the panels

library(patchwork)

TsfEcmSapVertical<- (TSFasvEcM1/TSFasvSap1)+plot_annotation(tag_levels = 'a');TsfEcmSapVertical
TsfEcmSapHorizontal<- (TSFasvEcM1|TSFasvSap1)+plot_annotation(tag_levels = 'a');TsfEcmSapHorizontal



# Export graphs-------------------------------------------------------------------------------------
dir.create("Analysis/MixedModels/PubGraphs")

pdf("Analysis/MixedModels/PubGraphs/TSF-vertical-EcM-Sap.pdf", height=11, width=7)
TsfEcmSapVertical
dev.off()

pdf("Analysis/MixedModels/PubGraphs/TSF-Hori-EcM-Sap.pdf", height=6, width=12)
TsfEcmSapHorizontal
dev.off()






