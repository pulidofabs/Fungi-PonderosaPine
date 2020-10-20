#SEPT 27, 2020
#TEST WHICH VARIABLES SIGNIFICANTLY INFLUENCE RICHNESS OF SOIL

#Reset R's Brain
rm(list=ls())

#* SET WORKING DIRECTORY.............................
#setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")

#* * LOAD LIBRARIES ********************************************************************----
library(lme4)#Linear generated mixed models analysis
library(lmerTest)#P-values for lme4 
library(ggplot2)
library(ggpubr)#For plots
library(ggfortify)#Data visualization


#*******************************************************************************************************************----
#* * LOAD DATA AND QUALITY CONTROL--------------------------------------------------------------------------------------
#*******************************************************************************************************************----
MetaRare<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobeNew.csv",na.strings = "N/A", header = TRUE) 


#Change reference level to Unburn 
MetaRare$Treatment<-as.factor(MetaRare$Treatment)
class(MetaRare$Treatment); levels(MetaRare$Treatment)
MetaRare$Treatment<- try(relevel(MetaRare$Treatment , "Unburned"))


#Scale data -----------------------------------------------------
MetaRare$Elevation<-scale(MetaRare$Elevation, scale = TRUE)
MetaRare$TC.TN<-scale(MetaRare$TC.TP, scale = TRUE)
MetaRare$TC.TP<-scale(MetaRare$TC.TN, scale = TRUE)  

head(MetaRare[1:2,])


#
#
#*************************************************************************************************************----
#* *BUILD MODEL (BACKWARD MODEL SELECTION)---------------------------------------------------------
#**************************************************************************************************************----

#TEST THE NULL MODEL ....................................................................................
attach(MetaRare)
#NULL model with no predictors 
Nullsp1<-lmer(S.obs~1+(1|Site/Plot/Subplot), data=MetaRare, REML = TRUE);summary(Nullsp1)#doees not work
Nullsp2<-lmer(S.obs~1+(1|Site/Plot), data=MetaRare, REML = TRUE);summary(Nullsp2)
Nullsp3<-lmer(S.obs~1+(1|Site), data=MetaRare, REML = TRUE);summary(Nullsp3)
Nullsp4<-lmer(S.obs~1+(1|Plot), data=MetaRare, REML = TRUE);summary(Nullsp4)

#Test the null model to see which one fits better .............
anova(Nullsp2,Nullsp3)#Null2sp 
anova(Nullsp2,Nullsp4)#Null4sp better



#Create the full model ....................................................................................
Fullsp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
               Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + 
               Treatment*TP +  Treatment*TN + Treatment*AvgDistanceTree +
               Treatment + TimeSinceFire + AvgDepthOM + TC.TN + TC.TP +
               TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
               Aspect + (1|Plot), MetaRare, REML = TRUE);summary(Fullsp) 
anova(Nullsp4,Fullsp)#Full model better 


#
#*************************************************************************************************************----
#* * BACKWARD MODEL SELECTION ---------------------------------------------------------
#**************************************************************************************************************----
#BEGIN BY REMOVING INTERACTION TERMS

#REMOVE Treatment*AvgDistanceTree +
Mod1sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
              Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + 
              Treatment*TP +  Treatment*TN + 
              Treatment + TimeSinceFire + AvgDepthOM + TC.TN + TC.TP +
              TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
              Aspect + (1|Plot), MetaRare, REML = TRUE);summary(Mod1sp) 
anova(Fullsp, Mod1sp)#Full model better, not significant, remove


#REMOVE  Treatment*TN +
Mod2sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
               Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + 
               Treatment*TP +   Treatment + TimeSinceFire +
               AvgDepthOM + TC.TN + TC.TP + TC + TP + TN + AvgDistanceTree +
               Slope + Elevation + Aspect + (1|Plot), MetaRare,
               REML = FALSE);summary(Mod2sp) 
anova(Fullsp, Mod2sp)#Mod 2 better, not sign, remove


#REMOVE Treatment*TP + 
Mod3sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
               Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + 
               Treatment + TimeSinceFire + AvgDepthOM + TC.TN + 
               TC.TP + TC + TP + TN + AvgDistanceTree + Slope + 
               Elevation + Aspect + (1|Plot), MetaRare,
               REML = FALSE);summary(Mod3sp) 
anova(Mod2sp, Mod3sp)#Mod3 better, not significant, remove



#REMOVE Treatment*TC + 
Mod4sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
              Treatment*TC.TN + Treatment*TC.TP +  Treatment +
              TimeSinceFire + AvgDepthOM + TC.TN +  TC.TP + 
              TC + TP + TN + AvgDistanceTree + Slope +  Elevation +
              Aspect + (1|Plot), MetaRare,
             REML = FALSE);summary(Mod4sp) 
anova(Mod3sp, Mod4sp)# Mod 4 better, not signficant, remove 



#REMOVE Treatment*TC.TP + 
Mod5sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
               Treatment*TC.TN +  Treatment +TimeSinceFire + 
               AvgDepthOM + TC.TN +  TC.TP +  TC + TP + TN + 
               AvgDistanceTree + Slope +  Elevation + Aspect +
              (1|Plot), MetaRare, REML = FALSE);summary(Mod5sp) 
anova(Mod4sp, Mod5sp)#Same AIC, not significant, remove


#REMOVE Treatment*TC.TN + 
Mod6sp<-lmer(S.obs ~ Treatment * TimeSinceFire + Treatment*AvgDepthOM +
                Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
                TC + TP + TN + AvgDistanceTree + Slope + Elevation +
                Aspect + (1|Plot), MetaRare, REML = FALSE);summary(Mod6sp) 
anova(Mod5sp, Mod6sp)# Mod 6 better, not signficant, remove 


#REMOVE Treatment*AvgDepthOM +
Mod7sp<-lmer(S.obs ~ Treatment * TimeSinceFire + 
               Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
               TC + TP + TN + AvgDistanceTree + Slope + Elevation +
               Aspect + (1|Plot), MetaRare, REML = FALSE);summary(Mod7sp) 
anova(Mod6sp, Mod7sp)# Mod 6 better, not siginficant,remove

#REMOVE Treatment * TimeSinceFire +
Mod8sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
               TC + TP + TN + AvgDistanceTree + Slope + Elevation +
               Aspect + (1|Plot), MetaRare, REML = FALSE);summary(Mod8sp) 
anova(Mod6sp, Mod8sp)#Mod 6 better, significant, do not remove......................

#REMOVE Aspect +
Mod9sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
               TC + TP + TN + AvgDistanceTree + Slope + Elevation +
               Treatment * TimeSinceFire +(1|Plot), 
             MetaRare, REML = FALSE);summary(Mod9sp) 
anova(Mod6sp, Mod9sp)#Mod 6 better, significant, do not remove ....................

#REMOVE Elevation +
Mod10sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
                TC + TP + TN + AvgDistanceTree + Slope + 
                Treatment * TimeSinceFire + Aspect + (1|Plot), 
              MetaRare, REML = FALSE);summary(Mod10sp) 
anova(Mod6sp, Mod10sp)#Mod 6 better, sign, do not remove ..........................


#Slope + 
Mod11sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
                TC + TP + TN + AvgDistanceTree +  Treatment * TimeSinceFire +
                Aspect + Elevation + (1|Plot),  MetaRare, 
                REML = FALSE);summary(Mod11sp) 
anova(Mod6sp, Mod11sp)#Mod 6 better, sign, do not remove ..........................

#AvgDistanceTree +  
Mod12sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + TC.TP + 
                TC + TP + TN + Treatment * TimeSinceFire +
                Aspect + Elevation + Slope + (1|Plot),  MetaRare, 
              REML = FALSE);summary(Mod12sp) 
anova(Mod6sp, Mod12sp)#SMod 12 better, not significant, remove


#TN + 
Mod13sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + 
                TC.TP + TC + TP + Treatment * TimeSinceFire + Aspect + 
                Elevation + Slope + (1|Plot),  MetaRare, 
                 REML = FALSE);summary(Mod13sp) 
anova(Mod12sp, Mod13sp)#Mod 13 better,not signficant, remove 

#TC +
Mod14sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + 
                TC.TP + TP + Treatment * TimeSinceFire + Aspect + 
                Elevation + Slope + (1|Plot),  MetaRare, 
              REML = FALSE);summary(Mod14sp) 
anova(Mod13sp, Mod4sp)#Mod 13 better, not significant, remove


#TP +
Mod15sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + 
                TC.TP + Treatment * TimeSinceFire + Aspect + 
                Elevation + Slope + (1|Plot),  MetaRare, 
              REML = FALSE);summary(Mod15sp) 
anova(Mod13sp, Mod15sp)#mOD 15 better, sign, do not remove ......................

#TC.TP + 
Mod16sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + TC.TN + 
               Treatment * TimeSinceFire + Aspect + Elevation +
               Slope + (1|Plot),  MetaRare, REML = FALSE);summary(Mod16sp) 
anova(Mod15sp, Mod16sp)#Same AIC, not significant, remove 

#TC.TN +
Mod17sp<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM +  
                Treatment * TimeSinceFire + Aspect + Elevation +
                Slope + (1|Plot),  MetaRare, REML = FALSE);summary(Mod17sp) 
anova(Mod16sp, Mod17sp)#Mod 17 better, not significant, remove

# AvgDepthOM +  
Mod18sp<-lmer(S.obs ~  Treatment +TimeSinceFire + 
                 Treatment * TimeSinceFire + Aspect + Elevation +
                 Slope + (1|Plot),  MetaRare, REML = FALSE);summary(Mod18sp) 
anova(Mod17sp, Mod18sp)#Mod 17 better, not significant, remove

# TimeSinceFire +
Mod19sp<-lmer(S.obs ~  Treatment + Treatment * TimeSinceFire +
                Aspect + Elevation + Slope + (1|Plot),  
               MetaRare, REML = FALSE);summary(Mod19sp) 
anova(Mod17sp, Mod19sp)#Mod 17 better, not significant, remove

# Treatment +
Mod20sp<-lmer(S.obs ~   Treatment * TimeSinceFire +
                Aspect + Elevation + Slope + (1|Plot),  
              MetaRare, REML = FALSE);summary(Mod20sp) 
anova(Mod17sp, Mod20sp)#Mod 20 better, not significant, remove


-----------------------------------------------------------------------------
anova(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp,
      Mod8sp, Mod9sp, Mod10sp, Mod11sp, Mod12sp,Mod13sp, 
      Mod14sp, Mod15sp, Mod16sp, Mod17sp, Mod18sp,
      Mod19sp, Mod20sp)


library(AICcmodavg)
models <- list(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp, Mod8sp,
               Mod9sp, Mod10sp,  Mod11sp,Mod12sp,Mod13sp, Mod14sp, Mod15sp,
               Mod16sp, Mod17sp, Mod18sp, Mod19sp,Mod20sp)



model.names <- c('Mod1sp','Mod2sp','Mod3sp','Mod4sp','Mod5sp','Mod6sp',
                 'Mod7sp','Mod8sp', 'Mod9sp','Mod10sp','Mod11sp','Mod12sp',
                 'Mod13sp','Mod14sp','Mod15sp', 'Mod16sp','Mod17sp','Mod18sp',
                 'Mod19sp', 'Mod20sp')

AICmodels<-aictab(cand.set = models, modnames = model.names);AICmodels

#Export model results........................................................................
write.csv(AICmodels, "Analysis/MixedModels/Ectomycorrhizal/LMER-SppRichness-Results.csv")

#BEST MODEL 
#-------- K ---AICc---Delta_AICc--AICcWt--Cum.Wt----LL
#Mod13sp 17  1865.16   0.00        0.13    0.13   -913.97

#
####
#*************************************************************************************************-----
# RUN FINAL MODEL -------------------------------------------------------------------------------------
#*************************************************************************************************-----
library(lmerTest)#Gives us P values 
library(MASS)
library(rcompanion)
library(MuMIn)

#FINAL MODEL .................................................................................
FinalModelsp<-lmer(S.obs ~   TC + TP + TN + AvgDistanceTree + Slope + Elevation + 
              Aspect + Treatment * TimeSinceFire + Treatment*AvgDepthOM + 
              AvgDepthOM + (1|Plot), MetaRare, REML = FALSE);summary(FinalModelsp)


#apply REML Mlmer()D TO FINAL MODEL 
MFinalSp<-lmer(S.obs ~   TC + TP + TN + AvgDistanceTree + Slope + Elevation + 
                Aspect + Treatment * TimeSinceFire + Treatment*AvgDepthOM + 
              AvgDepthOM + (1|Plot), MetaRare,  REML = TRUE);summary(MFinalSp)



#Get the R-square for the final model and the p-value for the variables
r.squaredGLMM(MFinalSp) 
summary(MFinalSp)

#EVALUATE MODEL 
Residuals<-resid(MFinalSp)


#******************************************************************************************----
#---------- CREATE PUBLICATION TABLES --------------------------------------------------------
#******************************************************************************************----
#Creat publication tables comparing best model #13 w lowest AIC and model 20, last model
library(sjPlot) # table functions

Mod13spFinal<-lmer(S.obs ~  Treatment +TimeSinceFire + AvgDepthOM + 
                TC.TN + TC.TP + TC + TP + Treatment * TimeSinceFire +
                Aspect +  Elevation + Slope + (1|Plot),  MetaRare,
                REML = TRUE);summary(Mod13spFinal)

Mod20spFinal<- lmer(S.obs ~ Treatment * TimeSinceFire +
                Aspect + Elevation + Slope + (1|Plot),  
                MetaRare, REML = TRUE);summary(Mod20spFinal)

Top2Mods<-tab_model(Mod13spFinal,Mod20spFinal, p.val = "kr", show.stat = TRUE);Top2Mods

#not sure how to export the table will copy and paste into excel















#
###
#************************************************************************************************************----
#------------------------   EVALUATE MODEL FIT ------------------------------------------------------------------
#************************************************************************************************************----
#This residual plot does not indicate any deviation from a linear form. It also shows relatively constant
# variance across the fitted range

plot(MFinalSp)


#Normality of residuals ..............................................
qqnorm(residuals(MFinalSp))


#Leverage plot........................................................
ggplot(data.frame(lev=hatvalues(MFinalSp),
                  pearson=residuals(MFinalSp,type="pearson")),
       aes(x=lev,y=pearson)) +geom_point() + theme_bw()


#prediction plot
newDat <- data.frame(x1=c(0.1,0.4),x2=c(-10,20))
predict(MFinalSp,newDat,re.form=NA,type="link")
predict(MFinalSp,newDat,re.form=NA,type="response")

newDat2 <- data.frame(x1=c(0.1,0.4),x2=c(-10,20),g1=c("B","B"))
predict(MFinalSp,newDat2,type="link")
predict(MFinalSp,newDat2,type="response")



#
##
#************************************************************************************************************----
#EXPORT MODEL RESULTS ------------------------------------------------------------------------------
#************************************************************************************************************----
dir.create(file.path("Analysis/MixedModels/Saprobes"), recursive = TRUE)
capture.output(summary(MFinalSp),file="Analysis/MixedModels/Saprobes/ASV-MixedModel-Results.csv")













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









