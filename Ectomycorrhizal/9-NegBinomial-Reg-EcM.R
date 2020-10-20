#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")


#LOAD LIBRARIES .................................................................
library(ggplot2)
library(tidyr)
library(dplyr)
library(ochRe)
library(MASS)
library(rcompanion)# for pseudo R


#LOAD DATA....................................................................... 
MetaRareEcM<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv",na.strings = "N/A", header = TRUE) 


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

levels(BurnedEcM$FireYear)


#Check the means and variance of the data to verify that neg binom is right
#variance is higher that mean= negative binomial 
df2 <- BurnedEcM %>%
  group_by(FireYear) %>%
  summarise(mcount = mean(S.obs),
            varcount = var(S.obs));df2


#RUN REGRESSION ---------------------------------------------------------------------------------------
library(MuMIn) #to get pseudo R
library(lme4)

isSingular(x, tol = 1e-4)

#Test if nested by Site, Plot, or subplot or combinations
ASVneg1 <- glmer.nb(S.obs ~ TimeSinceFire + (1|Site/Plot/Subplot),data = BurnedEcM);ASVneg1#does not converge
ASVneg2 <- glmer.nb(S.obs ~ TimeSinceFire + (1|Site/Plot),data = BurnedEcM);ASVneg2#does not converge
ASVneg3 <- glmer.nb(S.obs ~ TimeSinceFire + (1|Site),data = BurnedEcM);ASVneg3
ASVneg4 <- glmer.nb(S.obs ~ TimeSinceFire + (1|Plot),data = BurnedEcM);ASVneg4

anova(ASVneg1,ASVneg2)#ASVneg2 is better
anova(ASVneg2,ASVneg3)#ASVneg3 is better
anova(ASVneg3,ASVneg4)#ASVneg3 is better


# FINAL MODEL ----------------------- SITE AS RANDOM EFFECT ---------------------------------------
library(emmeans)

ASVnull <- glmer.nb(S.obs ~ 1 +( 1 | Site), data = BurnedEcM);summary(ASVnull)
anova(ASVneg, ASVnull)

#Get overall P value
ASVneg <- glmer.nb(S.obs ~ as.factor(TimeSinceFire) + (1|Site),data = BurnedEcM);summary(ASVneg)


#Redo as factor to do contrast
ASVneg1 <- glmer.nb(S.obs ~ as.factor(TimeSinceFire) + (1|Site),data = BurnedEcM);summary(ASVneg)
ASVemm <- emmeans(ASVneg1, "TimeSinceFire"); ASVemm
pairs(ASVemm)
pwpp(ASVemm)

Res<-summary(ASVneg);Res #0.00628 y = 2.29664X - 0.03492
AovRes<-anova(ASVneg); AovRes# F=7.4555

#CALCULATE R2 FOR THE MODEL -------------------------------------------------------
#R2M = marginal: for fixed effect  and  R2c= conditional for entire model 
R<-r.squaredGLMM(ASVneg);R





#PREDICTED GRAPHS------------------------------------------------------------------ 
library(ggeffects)
pr1 <- ggpredict(ASVneg, c("TimeSinceFire"))
plot(pr1)


#EXPORT RESULTS -----------------------------------------------------------------------------------------
dir.create(file.path("Analysis/MixedModels/Ectomycorrhizal/NegativeBinom"), recursive = TRUE)

capture.output(Res, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto.csv")
capture.output(AovRes, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto-ANOVA.csv")
capture.output(R, file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/Neg-Binom-Ecto-Rsquare.csv")
capture.output(summary(ASVneg), file="Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/ASV-NegBinom-ModelResults.csv")

##
##
#*********************************************************************************************************----
#--------------- MODEL VALIDATION ----------------------------------------------------------------------------
#*********************************************************************************************************----

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
library(ggpubr)

MetaRareEcM$FireYear<-as.factor(MetaRareEcM$TimeSinceFire)

TSFasvEcM<-ggplot(UnburnedEcM, aes(x=TimeSinceFire, y=S.obs))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE,aes(color=Treatment))+
  geom_point(aes(fill=factor(FireYear)), size=4, shape=21, stroke=0) +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
        sep = "~~~~")), formula = y~x, label.x = 1, label.y = 32, size=6)+
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
  labs(x = "Time Since Fire (Years)",y="Species Richness");TSFasvEcM


TSFasvEcM1<-TSFasvEcM + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  guides(color = guide_legend(override.aes = list(size=4))) ;TSFasvEcM1


  


#
#*********************************************************************************************************----
# --- EXPORT GRAPHS ------------------------------------------------------------------------------------------
#*********************************************************************************************************----


pdf("Analysis/MixedModels/Ectomycorrhizal/NegativeBinom/ASV-NB-TSF-EcM-Trt.pdf", height=6, width=7)
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

pdf("Analysis/MixedModels/PubGraphs/TSF-Trt-vertical-EcM-Sap.pdf", height=11, width=11)
TsfEcmSapVertical
dev.off()

pdf("Analysis/MixedModels/PubGraphs/TSF-Trt-Hori-EcM-Sap.pdf", height=8, width=14)
TsfEcmSapHorizontal
dev.off()


v



