#JUNE 14, 2020


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



#**************************************************************************************----
#* * LOAD DATA AND QUALITY CONTROL---------------------------------------------------------
#**************************************************************************************----

MetaRare<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.csv",
                    na.strings = "N/A", header = TRUE) 


#Change reference level to Unburn 
MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment<- try(relevel(MetaRare$Treatment , "Unburned"))


#Scale data -----------------------------------------------------
MetaRare$Elevation<-scale(MetaRare$Elevation, scale = TRUE)
MetaRare$TC.TN<-scale(MetaRare$TC.TP, scale = TRUE)
MetaRare$TC.TP<-scale(MetaRare$TC.TN, scale = TRUE)  

head(MetaRare[1:2,])

# SUBSET DATA ---------------------------------------------------
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

Null1sp<-lmer(S.obs~1+(1|Plot), data=MetaRare, REML = TRUE);summary(Null1sp)
anova(Nullsp,Null1sp)#Null1sp fits better


Fullsp<-lmer(S.obs ~ Treatment+Treatment*AvgDepthOM + AvgDepthOM+ Treatment*TC + TC + Treatment*TC.TN +
            TC.TN + Treatment*TP + TP + Treatment*TC.TP + TC.TP + Treatment*TN + TN + Slope + Elevation + 
            Aspect + AvgDistanceTree + (1|Plot), MetaRare, REML = TRUE);summary(Fullsp) 
anova(Fullsp,Null1sp)#Fullsp better 



#BACKWARD MODELS SELECTION------------------------------------------
#BEGIN BY REMOVING INTERACTION TERMS

#Remove Treatment*TC.TN
Mod1sp<-lmer(S.obs ~ TimeSinceFire+Treatment+Treatment*AvgDepthOM + AvgDepthOM+ 
          Treatment*TC + TC + TC.TN + Treatment*TP + TP + Treatment*TC.TP + TC.TP + 
          Treatment*TN + TN + Slope + Elevation + Aspect + AvgDistanceTree + (1|Plot), 
          MetaRare, REML = FALSE);summary(Mod1sp)
anova(Nullsp,Mod1sp)#Mod 1 better,  sign, add back in-----------------------------------------------------

#Remove Treatment*TC.TP 
Mod2sp<-lmer(S.obs ~ TimeSinceFire + Treatment + Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + Treatment*TC.TN + TC.TN +Treatment*TP + TP +
               TC.TP + Treatment*TN + TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod2sp)
anova(Mod1sp,Mod2sp)#Same AIC, not significant, remove

#Remove  Treatment*TN 
Mod3sp<-lmer(S.obs ~ TimeSinceFire + Treatment + Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + Treatment*TC.TN + TC.TN + Treatment*TP + TP +
               TC.TP + TN + Slope + Elevation + Aspect + 
               AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod3sp)
anova(Mod2sp,Mod3sp)#Mod 3 better, not significant, remove

#Remove  Treatment*Tp 
Mod4sp<-lmer(S.obs ~ TimeSinceFire + Treatment + Treatment*AvgDepthOM + AvgDepthOM+ 
               Treatment*TC + TC + Treatment*TC.TN + TC.TN + TP + TC.TP + TN +
               Slope + Elevation + Aspect + AvgDistanceTree + (1|Plot),
               MetaRare, REML = FALSE);summary(Mod4sp)
anova(Mod3sp,Mod4sp)#Mod 4 Better, Not significant, remove 


#Remove  Treatment*TC
Mod5sp<-lmer(S.obs ~ TimeSinceFire + Treatment + Treatment*AvgDepthOM + AvgDepthOM+ 
              TC + Treatment*TC.TN + TC.TN + TP + TC.TP + TN + Slope + Elevation +
              Aspect + AvgDistanceTree + (1|Plot), MetaRare, REML = FALSE);summary(Mod5sp)
anova(Mod4sp,Mod5sp)#Mod 5 better, not significant, remove


#Remove Treatment*AvgDepthOM
Mod6sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + Treatment*TC.TN +
               TC.TN + TP + TC.TP + TN + Slope + Elevation +Aspect + AvgDistanceTree + 
               (1|Plot), MetaRare, REML = FALSE);summary(Mod6sp)
anova(Mod5sp,Mod6sp)#Mod 5 better, not significant, remove

#Remove AvgDistanceTree 
Mod7sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + Treatment*TC.TN +
               TC.TN + TP + TC.TP + TN + Slope + Elevation +Aspect + (1|Plot), 
             MetaRare, REML = FALSE);summary(Mod7sp)
anova(Mod5sp,Mod7sp)#Mod 7 better, Not significant, remove


#remove Elevation
Mod8sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + Treatment*TC.TN +
               TC.TN + TP + TC.TP + TN + Slope + Aspect + (1|Plot), 
             MetaRare, REML = FALSE);summary(Mod8sp)
anova(Mod7sp,Mod8sp)#Mod 8 better, Not significant, remove


#remove Aspect 
Mod9sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + Treatment*TC.TN +
               TC.TN + TP + TC.TP + TN + Slope +  (1|Plot), 
             MetaRare, REML = FALSE);summary(Mod9sp)
anova(Mod8sp,Mod9sp)#Mod 8 better, significant, add back in ----------------------------------

#remove Slope
Mod10sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + 
              Treatment*TC.TN + TC.TN + TP + TC.TP + TN + Aspect + (1|Plot),
              MetaRare, REML = FALSE);summary(Mod10sp)
anova(Mod8sp,Mod10sp)#Mod 10 better, not significant remove

#remove TN
Mod11sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + 
                Treatment*TC.TN + TC.TN + TP + TC.TP +  Aspect + (1|Plot),
              MetaRare, REML = FALSE);summary(Mod11sp)
anova(Mod10sp,Mod11sp)#Mod 11 better, not sign, remove


#Remove Tp
Mod12sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM+ TC + 
                Treatment*TC.TN + TC.TN +  TC.TP +  Aspect + (1|Plot),
              MetaRare, REML = FALSE);summary(Mod12sp)
anova(Mod11sp,Mod12sp)#Mod 12 better, not significant,  remove

#remove TC
Mod13sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM +  
                Treatment*TC.TN + TC.TN +  TC.TP +  Aspect + (1|Plot),
              MetaRare, REML = FALSE);summary(Mod13sp)
anova(Mod12sp,Mod13sp)#MOD 12 BETTER, SIGN, add back in -------------------------------------------------

#remove TC.TP  
Mod14sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM +  
                Treatment*TC.TN + TC.TN + Aspect + TC +
                (1|Plot),  MetaRare, REML = FALSE);summary(Mod14sp)
anova(Mod12sp,Mod14sp)#Same AIC, not significant, remvoe


#REMOVE TC.TN
Mod15sp<-lmer(S.obs ~ TimeSinceFire + Treatment + AvgDepthOM +  
                Treatment*TC.TN + Aspect + TC + (1|Plot),  
              MetaRare, REML = FALSE);summary(Mod15sp)
anova(Mod14sp,Mod15sp)#Same AIC, not significant, remove


#remove AvgDepthOM+ 
Mod16sp<-lmer(S.obs ~ TimeSinceFire + Treatment +  Treatment*TC.TN + 
                Aspect + TC + (1|Plot), MetaRare, REML = FALSE);summary(Mod16sp)
anova(Mod15sp,Mod16sp)#Mod 16 better, not sign, remove


#remove TREATMENT
Mod17sp<-lmer(S.obs ~ TimeSinceFire +  Treatment*TC.TN + Aspect + TC +
                (1|Plot), MetaRare, REML = FALSE);summary(Mod17sp)
anova(Mod16sp,Mod17sp)#Same AIC, not significant, remove


#Remove TimeSinceFire
Mod18sp<-lmer(S.obs ~ Treatment*TC.TN + Aspect + TC + (1|Plot),
               MetaRare, REML = FALSE);summary(Mod18sp)
anova(Mod17sp,Mod18sp)#Mod 17 better, sign, add back in ---------------------------------





#FINAL MODEL--MOD 17
lmer(S.obs ~ TimeSinceFire +  Treatment*TC.TN + Aspect + TC +(1|Plot), MetaRare, REML = FALSE)



anova(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp, Mod8sp, Mod9sp, Mod10sp, Mod11sp,
      Mod12sp,Mod13sp, Mod14sp, Mod15sp, Mod16sp, Mod17sp, Mod18sp)


#*************************************************************************************************-----
# RUN FINAL MODEL -------------------------------------------------------------------------------------
#*************************************************************************************************-----
library(lmerTest)#Gives us P values 
library(MASS)
library(rcompanion)

#AIC -39.---- Trtm, Elev, Trt*TN sign p<0.05
FinalModelsp<-lmer(S.obs ~ TimeSinceFire + TC +  Treatment*TC.TN + Aspect + (1|Plot), 
                   MetaRare, REML = FALSE);summary(FinalModelsp)

#apply REML Mlmer()D TO FINAL MODEL 
MFinalSp<-lmer(S.obs ~ TimeSinceFire + TC +  Treatment*TC.TN + Aspect + 
                 (1|Plot), MetaRare, REML = TRUE);summary(MFinalSp)



MuMIn::r.squaredGLMM(lme4::lmer(S.obs ~ TimeSinceFire + TC +  Treatment*TC.TN + Aspect + 
                (1|Plot), MetaRare, REML = TRUE), RE.keep=TRUE)




#EVALUATE MODEL 
Residuals<-resid(MFinalSp)


#************************************************************************************************************----
#EXPORT MODEL RESULTS ------------------------------------------------------------------------------
#************************************************************************************************************----
dir.create("Analysis/MixedModels/Saprobes")
capture.output(summary(MFinalSp),file="Analysis/MixedModels/Saprobes/ASV-MixedModel-Results.csv")








#MAKE SURE EVERYTHING IS OKAY W THE MODEL-------------------------------------------------
E<-resid(MFinalSp)#extract residuals from the fitted model
Fit<-fitted(MFinalSp)#extract the fitted values
op<-par(mfrow=c(1,2))
plot(x=Fit, y=E, xlab="Fitted Values", 
     ylab="Residuals", main="Residuals vs. Fitted Values")
identity(Fit,E)#does not work
hist(E,nclass = 15)
par(op)






#
#
#*************************************************************************************************************----
#* *BUILD MODEL (BACKWARD MODEL SELECTION) BURNED---------------------------------------------------------
#**************************************************************************************************************----

#NULL model with no predictors 
Nullsp<-lmer(S.obs~1+(1|Site/Plot), data=Burned, REML = TRUE);summary(Nullsp)

Null1sp<-lmer(S.obs~1+(1|Site), data=Burned, REML = TRUE);summary(Null1sp)
anova(Nullsp,Null1sp)#Null1sp fits better


Fullsp<-lmer(S.obs ~ TimeSinceFire+AvgDepthOM+  TC + TC.TN +  TP + TC.TP + TN + Slope + 
               Elevation + Aspect + AvgDistanceTree + (1|Plot), Burned, REML = TRUE);summary(Fullsp) 
anova(Fullsp,Null1sp)#Fullsp better 



#BACKWARD MODELS SELECTION------------------------------------------
#BEGIN BY REMOVING INTERACTION TERMS

#AvgDepthOM
Mod1sp<-lmer(S.obs ~ TimeSinceFire + TC + TN +TP + TC.TN + TC.TP + 
               Slope + Elevation + Aspect + AvgDistanceTree + (1|Plot),
             Burned, REML = TRUE);summary(Mod1sp)
anova(Fullsp,Mod1sp)#Full Mod better, significant, do not remove



#Aspect
Mod2sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN +  TP +
               TC.TP + TN + Slope + Elevation + AvgDistanceTree +
               (1|Plot), Burned, REML = TRUE);summary(Mod2sp)
anova(Fullsp,Mod2sp)#Mod 2 better, not significant, keep off

# Elevation
Mod3sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN + TP +
               TC.TP + TN + Slope + AvgDistanceTree +
               (1|Plot), Burned, REML = TRUE);summary(Mod3sp)
anova(Mod2sp,Mod3sp)#Mod 2 better, Not significant, remove


#Slope 
Mod4sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN + TP +
               TC.TP + TN + AvgDistanceTree +
               (1|Plot), Burned, REML = TRUE);summary(Mod4sp)
anova(Mod2sp,Mod4sp)#Mod 4 better, Not significant, remove



#AvgDepthOM
Mod5sp<-lmer(S.obs ~ TimeSinceFire +  TC + TC.TN + TP +
               TC.TP + TN + AvgDistanceTree +
               (1|Plot), Burned, REML = TRUE);summary(Mod5sp)
anova(Mod4sp,Mod5sp)#Mod 4 better, significant, Add back in 


#AvgDistTree
Mod6sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN + TP +
               TC.TP + TN +  (1|Plot), Burned, REML = TRUE);summary(Mod6sp)
anova(Mod4sp,Mod6sp)#Mod 6 better, not significant remove

#TN
Mod7sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN + TP +
               TC.TP +  (1|Plot), Burned, REML = TRUE);summary(Mod7sp)
anova(Mod6sp,Mod7sp)#Mod 7 better, not significant remove


#TCTP
Mod8sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN + TP +
                (1|Plot), Burned, REML = TRUE);summary(Mod8sp)
anova(Mod7sp,Mod8sp)#Mod 8 better, not significant remove


#rEMOVE tp
Mod9sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC + TC.TN +
               (1|Plot), Burned, REML = TRUE);summary(Mod9sp)
anova(Mod8sp,Mod9sp)#Mod9 better, not significant remove

#rEMOVE TC.TN
Mod10sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+  TC +
               (1|Plot), Burned, REML = TRUE);summary(Mod10sp)
anova(Mod9sp,Mod10sp)#Mod10 better, not significant remove

#rEMOVE TC
Mod11sp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+ +
                (1|Plot), Burned, REML = TRUE);summary(Mod11sp)
anova(Mod10sp,Mod11sp)#Mod11 better, not significant remove


#TimeSinceFire
Mod12sp<-lmer(S.obs ~ AvgDepthOM + 
                (1|Plot), Burned, REML = TRUE);summary(Mod12sp)
anova(Mod11sp,Mod12sp)#Mod12 Significant, put back in 




anova(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp, Mod8sp, Mod9sp, Mod10sp, Mod11sp,
      Mod12sp)


#*************************************************************************************************-----
# RUN FINAL MODEL -------------------------------------------------------------------------------------
#*************************************************************************************************-----
library(lmerTest)#Gives us P values 

#AIC -39.---- Trtm, Elev, Trt*TN sign p<0.05
FinalModelsp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+ +
                     (1|Plot), Burned, REML = FALSE);summary(FinalModelsp)

#apply REML Mlmer()D TO FINAL MODEL 
MFinalSp<-lmer(S.obs ~ TimeSinceFire + AvgDepthOM+ +
                 (1|Plot), Burned, REML = TRUE);summary(MFinalSp)

summary(MFinalSp)



#EVALUATE MODEL 
Residuals<-resid(MFinalSp)


#************************************************************************************************************----
#EXPORT MODEL RESULTS ------------------------------------------------------------------------------
#************************************************************************************************************----
dir.create("Analysis/MixedModels/Ectomycorrhizal")
capture.output(summary(MFinalSp),file="Analysis/MixedModels/Ectomycorrhizal/ASV-BURNED-MixedModel-Results.csv")








#MAKE SURE EVERYTHING IS OKAY W THE MODEL-------------------------------------------------
E<-resid(MFinalSp)#extract residuals from the fitted model
Fit<-fitted(MFinalSp)#extract the fitted values
op<-par(mfrow=c(1,2))
plot(x=Fit, y=E, xlab="Fitted Values", 
     ylab="Residuals", main="Residuals vs. Fitted Values")
identity(Fit,E)#does not work
hist(E,nclass = 15)
par(op)









