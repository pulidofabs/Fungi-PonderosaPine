#TEST WHICH VARIABLES SIGNIFICANTLY INFLUENCE RICHNESS OF SOIL

#Reset R's Brain
rm(list=ls())

#* SET WORKING DIRECTORY.............................
setwd("C:/Users/fabipc/Dropbox/PIPO")

#* * LOAD LIBRARIES ********************************************************************----
library(lme4)#Linear generated mixed models analysis
library(lmerTest)#P-values for lme4 
library(ggplot2)
library(ggpubr)#For plots
library(ggfortify)#Data visualization
library(rcompanion)#pseudoR for regression



#* * LOAD DATA AND QUALITY CONTROL---------------------------------------------------------
MetaRare<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv",
                    na.strings = "N/A", header = TRUE) 

#***************************************************************************************----
# QUALITY CONTROL --------------------------------------------------------------------------
#***************************************************************************************---- 

#Change reference level to Unburn ................................................
MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment<- try(relevel(MetaRare$Treatment , "Unburned"))


#Scale data ......................................................................
MetaRare$Elevation<-scale(MetaRare$Elevation, scale = TRUE)
MetaRare$TC.TN<-scale(MetaRare$TC.TP, scale = TRUE)
MetaRare$TC.TP<-scale(MetaRare$TC.TN, scale = TRUE)  



# SUBSET DATA --------------------------------------------------------------------
# # * * UNBURNED PLOTS ------------------------------
attach(MetaRare)
Unburned<-MetaRare[which(Treatment== "Unburned"), ]
head(Unburned);dim(Unburned)#105x33

# * * BURNED PLOTS ------------------------------
Burned<-MetaRare[which(Treatment== "Burned"), ]
head(Burned);dim(Burned)#87x33
detach(MetaRare)


#
#
#*************************************************************************************************************----
#* *BUILD MODEL (BACKWARD MODEL SELECTION)---------------------------------------------------------
#**************************************************************************************************************----

#NULL model with no predictors 
Nullsp<-lmer(S.obs~1+(1|Site/Plot), data=MetaRare, REML = TRUE);summary(Nullsp)

Null1sp<-lmer(S.obs~1+(1|Site), data=MetaRare, REML = TRUE);summary(Null1sp)
anova(Nullsp,Null1sp)#Null1sp fits better


Fullsp<-lmer(S.obs ~ Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + Treatment*TC.TN +
            TC.TN + Treatment*TP + TP + Treatment*TC.TP + TC.TP + Treatment*TN + TN + Slope + Elevation + 
            Aspect + AvgDistanceTree + (1|Plot), MetaRare, REML = TRUE);summary(Fullsp) 
anova(Fullsp,Null1sp)#Fullsp better 



#BACKWARD MODELS SELECTION------------------------------------------
#BEGIN BY REMOVING INTERACTION TERMS

#Remove Treatment*TC.TN
Mod1sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + TC.TN +
               Treatment*TP + TP + Treatment*TC.TP + TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod1sp)
anova(Fullsp,Mod1sp)#Mod 1 better, but variable sign, add back in

#Remove Treatment*TC.TP 
Mod2sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + TC.TN +
               Treatment*TP + TP +  TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod2sp)
anova(Mod1sp,Mod2sp)#Mod 2 better, Not significant, remove


#Remove  Treatment*TN 
Mod3sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + TC.TN +
               Treatment*TP + TP +  TC.TP +  TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod3sp)
anova(Mod2sp,Mod3sp)#Mod 2, significant, do not remove

#Remove  Treatment*Tp add back in Trt * TN
Mod4sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + TC.TN +
               TP +  TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
              AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod4sp)
anova(Mod2sp,Mod4sp)#Mod 4 better, not sign, remove


#Remove  Treatment*TC
Mod5sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM + TC + TC.TN +
               TP +  TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod5sp)
anova(Mod4sp,Mod5sp)#Mod 4 better, Treatment and TC significant, do not remove


#Remove Treatment*AvgDepthOM and add back in Trt * TC
Mod6sp<-lmer(S.obs ~ TimeSinceFire+Treatment+AvgDepthOM+ Treatment*TC + TC + TC.TN +
                TP +  TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
                AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod6sp)
anova(Mod4sp,Mod6sp)#Mod 4 better, significant, do not remove


#remove Elevation
Mod7sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + TC.TN + TP +  TC.TP + Treatment*TN + TN + Slope + Aspect + 
              AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod7sp)
anova(Mod4sp,Mod7sp)#Mod 7 better, Not significant, remove


#Remove AvgDistanceTree 
Mod8sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + TC.TN + TP +  TC.TP + Treatment*TN + TN + Slope + Aspect + 
               (1|Plot), MetaRare, REML = FALSE);summary(Mod8sp)
anova(Mod7sp,Mod8sp)#Mod 8 better, not sign, remove


#remove Aspect 
Mod9sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + TC.TN + TP +  TC.TP + Treatment*TN + TN + Slope +
               (1|Plot), MetaRare, REML = FALSE);summary(Mod9sp)
anova(Mod8sp,Mod9sp)#Mod 9 better, Not significant, remove

#remove Slope
Mod10sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                Treatment*TC + TC + TC.TN + TP +  TC.TP + Treatment*TN + TN +
                (1|Plot), MetaRare, REML = FALSE);summary(Mod10sp)
anova(Mod9sp,Mod10sp)#Mod 9 better, not significant, remove


#Remove TN
Mod11sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                Treatment*TC + TC + TC.TN + TP +  TC.TP + Treatment*TN + 
                (1|Plot), MetaRare, REML = FALSE);summary(Mod11sp)
anova(Mod9sp,Mod11sp)#Mod 9 bette, not sign, remove

#Remove Tp
Mod12sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                Treatment*TC + TC + TC.TN + TC.TP + Treatment*TN + 
                (1|Plot), MetaRare, REML = FALSE);summary(Mod12sp)
anova(Mod9sp,Mod12sp)#Mod 12 better, Not significant remove

#remove TC
Mod13sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                 Treatment*TC + TC.TN + TC.TP + Treatment*TN + 
                 (1|Plot), MetaRare, REML = FALSE);summary(Mod13sp)
anova(Mod12sp,Mod13sp)#Same AIC, Not significant remove

#remove TC.TP  
Mod14sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                Treatment*TC + TC.TN + Treatment*TN + 
                (1|Plot), MetaRare, REML = FALSE);summary(Mod14sp)
anova(Mod13sp,Mod14sp)#Same AIC, not significant remove


#REMOVE TC.tn
Mod15sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
                Treatment*TC + Treatment*TN + 
                (1|Plot), MetaRare, REML = FALSE);summary(Mod15sp)
anova(Mod14sp,Mod15sp)#Mod 14 better, Not significant, remove

#remove AvgDepthOM+ 
Mod16sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM +
                Treatment*TC + Treatment*TN + (1|Plot), MetaRare, REML = FALSE);summary(Mod16sp)
anova(Mod14sp,Mod16sp)#Mod 14 better, not sign remove


#remove Treatment
Mod17sp<-lmer(S.obs ~ TimeSinceFire+Treatment*AvgDepthOM +
                Treatment*TC + Treatment*TN + (1|Plot), MetaRare, REML = FALSE);summary(Mod17sp)
anova(Mod14sp,Mod17sp)#Mod 14 better, not sign, remove

#Remove TSF 
Mod18sp<-lmer(S.obs ~ Treatment*AvgDepthOM +Treatment*TC + Treatment*TN + (1|Plot), 
              MetaRare, REML = FALSE);summary(Mod18sp)
anova(Mod14sp,Mod18sp)#Mod 14 better, TSF Significant, do not remove



#FINAL MODEL--model 14 better
anova(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp, Mod8sp, Mod9sp, Mod10sp, Mod11sp,
      Mod12sp,Mod13sp, Mod14sp, Mod15sp, Mod16sp, Mod17sp, Mod18sp)


#*************************************************************************************************-----
# RUN FINAL MODEL -------------------------------------------------------------------------------------
#*************************************************************************************************-----
library(lmerTest)#Gives us P values 
library(MuMIn)

#AIC -39.---- Trtm, Elev, Trt*TN sign p<0.05
FinalModelsp<-lmer(S.obs ~ TimeSinceFire+Treatment*AvgDepthOM + Treatment*TC + 
                    Treatment*TN + (1|Plot), MetaRare, REML = FALSE);summary(FinalModelsp)

#apply REML Mlmer()D TO FINAL MODEL 
MFinalSp<-lmer(S.obs ~ TimeSinceFire+Treatment*AvgDepthOM + Treatment*TC + 
                 Treatment*TN + (1|Plot), MetaRare,  REML = TRUE);summary(MFinalSp)


MuMIn::r.squaredGLMM(lme4::lmer(S.obs ~ TimeSinceFire+Treatment*AvgDepthOM + Treatment*TC + 
                                  Treatment*TN + (1|Plot), MetaRare,  REML = TRUE), RE.keep=TRUE)
summary(MFinalSp)


#EVALUATE MODEL-----------------------------
Residuals<-resid(MFinalSp)




#************************************************************************************************************----
#EXPORT MODEL RESULTS ------------------------------------------------------------------------------
#************************************************************************************************************----
dir.create("Analysis/MixedModels/Ectomycorrhizal")
capture.output(summary(MFinalSp),file="Analysis/MixedModels/Ectomycorrhizal/ASV-MixedModel-Results.csv")









