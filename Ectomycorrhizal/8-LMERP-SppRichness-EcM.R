#TEST WHICH VARIABLES SIGNIFICANTLY INFLUENCE RICHNESS OF SOIL
#Sept 27

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
library(rcompanion)#pseudoR for regression



#* * LOAD DATA AND QUALITY CONTROL---------------------------------------------------------
MetaRare<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv",
                    na.strings = "NA", header = TRUE) 

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
attach(MetaRare)
#NULL model with no predictors 
Nullsp1<-lmer(S.obs~1+(1|Site/Plot/Subplot), data=MetaRare, REML = TRUE);summary(Nullsp1)
Nullsp2<-lmer(S.obs~1+(1|Site/Plot), data=MetaRare, REML = TRUE);summary(Nullsp2)
Nullsp3<-lmer(S.obs~1+(1|Site), data=MetaRare, REML = TRUE);summary(Nullsp3)
Nullsp4<-lmer(S.obs~1+(1|Plot), data=MetaRare, REML = TRUE);summary(Nullsp4)

#Test the null model to see which one fits better .............
anova(Nullsp1,Nullsp2)#Null3sp better
anova(Nullsp3,Nullsp4)#Null3sp better



#Create full model to compare to null model .......................................................................
Fullsp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment*AvgDepthOM + Treatment*TC.TN +
               Treatment*TC.TP + Treatment*TC + Treatment*TP +  Treatment*TN + 
               Treatment*AvgDistanceTree + Treatment + TimeSinceFire + AvgDepthOM + 
               TC.TN + TC.TP + TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
               Aspect + (1|Site), MetaRare, REML = TRUE);summary(Fullsp) 
anova(Fullsp,Nullsp1)#Same AIC, proceed w full model 


#
#*************************************************************************************************************----
#* * BACKWARD MODEL SELECTION ---------------------------------------------------------
#**************************************************************************************************************----
#BEGIN BY REMOVING INTERACTION TERMS


#Remove Treatment*AvgDistanceTree
Mod1sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
          Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + Treatment*TP + 
          Treatment*TN + Treatment + TimeSinceFire + AvgDepthOM +  TC.TN + 
          TC.TP + TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
          Aspect + (1|Site), MetaRare, REML = FALSE);summary(Mod1sp) 
anova(Fullsp, Mod1sp)#Mod 1 better, not significant, remove

#Remove Treatment*TN + 
Mod2sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
        Treatment*TC.TN + Treatment*TC.TP + Treatment*TC + Treatment*TP +  
        Treatment + TimeSinceFire + AvgDepthOM + TC.TN + TC.TP +  TC +
        TP + TN +  AvgDistanceTree + Slope + Elevation + Aspect + 
        (1|Site), MetaRare, REML = FALSE);summary(Mod2sp) 
anova(Mod1sp, Mod2sp)#Mod 1 better, significant, do not remove ......................

#Remove Treatment*TP + 
Mod3sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
         Treatment*TC.TN + Treatment*TC.TP + Treatment*TC +  
         Treatment + TimeSinceFire + AvgDepthOM + TC.TN + TC.TP +  TC +
         TP + TN +  AvgDistanceTree + Slope + Elevation + Aspect + 
         Treatment*TN + (1|Site), MetaRare, REML = FALSE);summary(Mod3sp) 
anova(Mod1sp, Mod3sp)#Mod 3 better, not significant, remove ........................


#Remove Treatment*TC + 
Mod4sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
          Treatment*TC.TN + Treatment*TC.TP + Treatment + TimeSinceFire +
          AvgDepthOM + TC.TN + TC.TP +  TC + TP + TN +  AvgDistanceTree + 
          Slope + Elevation + Aspect + Treatment*TN + (1|Site), MetaRare,
          REML = FALSE);summary(Mod4sp) 
anova(Mod3sp, Mod4sp)#Mod 3 better, significant, do not remove......................


#Remove Treatment*TP + 
Mod5sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
          Treatment*TC.TN + Treatment*TC.TP + Treatment + TimeSinceFire +
          AvgDepthOM + TC.TN + TC.TP +  TC + TP + TN +  AvgDistanceTree + 
          Slope + Elevation + Aspect + Treatment*TN + Treatment*TC + 
          (1|Site), MetaRare,REML = FALSE);summary(Mod5sp)
anova(Mod3sp, Mod5sp)#Same AIC, not significant, remove



#Remove Treatment*TC.TP +
Mod6sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
        Treatment*TC.TN +  Treatment + TimeSinceFire +
        AvgDepthOM + TC.TN + TC.TP +  TC + TP + TN +  AvgDistanceTree + 
        Slope + Elevation + Aspect + Treatment*TN + Treatment*TC + 
        (1|Site), MetaRare,REML = FALSE);summary(Mod6sp) 
anova(Mod5sp, Mod6sp)#Same AIC, not significant, remove


#Remove Treatment*TC.TN +
Mod7sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment * AvgDepthOM +
          Treatment + TimeSinceFire + AvgDepthOM + TC.TN + TC.TP +
          TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
          Aspect +Treatment*TN + Treatment*TC + (1|Site), 
          MetaRare,REML = FALSE);summary(Mod7sp) 
anova(Mod6sp, Mod7sp)#Mod 7 sligtly better, not significant, remove

#Remove Treatment * AvgDepthOM +
Mod8sp<-lmer(S.obs~Treatment * TimeSinceFire + Treatment + TimeSinceFire +
          AvgDepthOM + TC.TN + TC.TP + TC + TP + TN +  AvgDistanceTree +
          Slope + Elevation + Aspect +Treatment*TN + Treatment*TC + 
          (1|Site), MetaRare,REML = FALSE);summary(Mod8sp) 
anova(Mod7sp, Mod8sp)#Mod 7 better, significant, do not remove .....................

#Remove Treatment * TimeSinceFire +
Mod9sp<-lmer(S.obs~ Treatment + TimeSinceFire +
         AvgDepthOM + TC.TN + TC.TP + TC + TP + TN +  AvgDistanceTree +
         Slope + Elevation + Aspect +Treatment*TN + Treatment*TC + 
         Treatment * AvgDepthOM +(1|Site), MetaRare,REML = FALSE);summary(Mod9sp) 
anova(Mod7sp, Mod9sp)#Mod 9 better, not significant, remove 

#Remove Aspect +
Mod10sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM + TC.TN +
          TC.TP + TC + TP + TN +  AvgDistanceTree + Slope + Elevation + 
          Treatment*TN + Treatment*TC + Treatment * AvgDepthOM +
          (1|Site), MetaRare,REML = FALSE);summary(Mod10sp) 
anova(Mod9sp, Mod10sp)#Mod 10 better, not significant, remove 

# Remove Elevation +
Mod11sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM + TC.TN +
                TC.TP + TC + TP + TN +  AvgDistanceTree + Slope +  
                Treatment*TN + Treatment*TC + Treatment * AvgDepthOM +
                (1|Site), MetaRare,REML = FALSE);summary(Mod11sp) 
anova(Mod10sp, Mod11sp)#Mod 11 better, not significant, remove 


# Remove Slope +
Mod12sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM + TC.TN +
                TC.TP + TC + TP + TN +  AvgDistanceTree +   
                Treatment*TN + Treatment*TC + Treatment * AvgDepthOM +
                (1|Site), MetaRare,REML = FALSE);summary(Mod12sp) 
anova(Mod11sp, Mod12sp)#Mod 11 slightly better (.2) not signficant, remove

# Remove  AvgDistanceTree +
Mod13sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN + TC.TP + TC + TP + TN + Treatment*TN +    
              Treatment*TC + Treatment * AvgDepthOM +(1|Site),
              MetaRare,REML = FALSE);summary(Mod13sp) 
anova(Mod11sp, Mod13sp)#Model 13 better, not significant, remove 


# Remove TN
Mod14sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN + TC.TP + TC + TP + Treatment*TN +    
                Treatment*TC + Treatment * AvgDepthOM +(1|Site),
              MetaRare,REML = FALSE);summary(Mod14sp) 
anova(Mod13sp, Mod14sp)#Same AIC, not significant, remove


# Remove TP + 
Mod15sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN + TC.TP + TC + Treatment*TN +    
                Treatment*TC + Treatment * AvgDepthOM +(1|Site),
              MetaRare,REML = FALSE);summary(Mod15sp) 
anova(Mod14sp, Mod15sp)#Mod 15 better, not significant, remove


# Remove TC + 
Mod16sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN + TC.TP + Treatment*TN +Treatment*TC + 
                Treatment * AvgDepthOM +(1|Site),MetaRare,
                REML = FALSE);summary(Mod16sp) 
anova(Mod15sp, Mod16sp)#Same AIC, not significant, remove


# Remove  TC.TP +
Mod17sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN +  Treatment*TN +Treatment*TC + 
                Treatment * AvgDepthOM +(1|Site),MetaRare,
              REML = FALSE);summary(Mod17sp) 
anova(Mod16sp, Mod17sp)#Same AIC, not significant, remove


# Remove  TC.TN + 
Mod18sp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
          Treatment*TN +Treatment*TC + Treatment * AvgDepthOM +
          (1|Site),MetaRare, REML = FALSE);summary(Mod18sp) 
anova(Mod17sp, Mod18sp)#Mod 17 better, not significant, remove

# Remove AvgDepthOM +
Mod19sp<-lmer(S.obs~ Treatment + TimeSinceFire +  
                Treatment*TN +Treatment*TC + Treatment * AvgDepthOM +
                (1|Site),MetaRare, REML = FALSE);summary(Mod19sp) 
anova(Mod17sp, Mod19sp)#Mod 17 better, not sigificant, remove 


# Remove TimeSinceFire +  
Mod20sp<-lmer(S.obs~ Treatment + Treatment*TN +Treatment*TC + 
          Treatment * AvgDepthOM + (1|Site),MetaRare, 
          REML = FALSE);summary(Mod20sp) 
anova(Mod17sp, Mod20sp)#Mod 17 better, significant, do not remove.................

#Remove Treatment +
Mod21sp<-lmer(S.obs ~  Treatment*TN +Treatment*TC + 
              Treatment * AvgDepthOM + TimeSinceFire +(1|Site),
              MetaRare, REML = FALSE);summary(Mod21sp) 
anova(Mod17sp, Mod21sp)#Mod 17 better, not significant, remove



#
###
#********************************************************************************************************------
# --------------- COMPARE MODELS TO SELECT FINAL MODEL --------------------------------------------------------
#********************************************************************************************************------
#COMPARE ALL THE MODELS TO SEE WHICH ONE IS BETTER 
anova(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp,
      Mod8sp, Mod9sp, Mod10sp, Mod11sp,Mod12sp,Mod13sp,
       Mod14sp, Mod15sp, Mod16sp, Mod17sp, Mod18sp, Mod19sp,
      Mod20sp, Mod21sp)

library(AICcmodavg)
models <- list(Mod1sp,Mod2sp,Mod3sp,Mod4sp, Mod5sp, Mod6sp, Mod7sp, Mod8sp,
               Mod9sp, Mod10sp,  Mod11sp,Mod12sp,Mod13sp, Mod14sp, Mod15sp,
               Mod16sp, Mod17sp, Mod18sp, Mod19sp,Mod20sp,Mod21sp)



model.names <- c('Mod1sp','Mod2sp','Mod3sp','Mod4sp','Mod5sp','Mod6sp',
                 'Mod7sp','Mod8sp', 'Mod9sp','Mod10sp','Mod11sp','Mod12sp',
                 'Mod13sp','Mod14sp','Mod15sp', 'Mod16sp','Mod17sp','Mod18sp',
                 'Mod19sp', 'Mod20sp','Mod21sp')

AICmodels<-aictab(cand.set = models, modnames = model.names);AICmodels

#Export model results........................................................................
write.csv(AICmodels, "Analysis/MixedModels/Ectomycorrhizal/LMER-SppRichness-Results.csv")


#
##
#*************************************************************************************************-----
# RUN FINAL MODEL -------------------------------------------------------------------------------------
#*************************************************************************************************-----
library(lmerTest)#Gives us P values 
library(MuMIn)

#Best models were 15, 16 and 17 based on AIC and AICccweight, will report the 
#results of the 3 models on table but will focus on model 17 as the final model


FinalModelsp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN +Treatment*TN +Treatment*TC + Treatment * AvgDepthOM +
                (1|Site),MetaRare,REML = FALSE);summary(FinalModelsp)


#apply REML Mlmer()D TO FINAL MODEL 
MFinalSp<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM + TC.TN +
                 Treatment*TN +Treatment*TC + Treatment * AvgDepthOM +
                 (1|Site),MetaRare,REML = TRUE);summary(MFinalSp)


#Get the R-square for the final model and the p-value for the variables
r.squaredGLMM(MFinalSp) 
summary(MFinalSp)


#----
#******************************************************************************************----
#---------- CREATE PUBLICATION TABLES --------------------------------------------------------
#******************************************************************************************----
#Creat publication tables
library(sjPlot) # table functions

#Put the results of the top 3 models together
Mod17spFinal<-lmer(S.obs~ Treatment + TimeSinceFire +  AvgDepthOM +
                TC.TN +  Treatment*TN +Treatment*TC + 
                Treatment * AvgDepthOM + (1|Site),MetaRare,
                REML = TRUE);summary(Mod17spFinal) 

Mod21spFinal<-lmer(S.obs ~  Treatment*TN +Treatment*TC + 
              Treatment * AvgDepthOM + TimeSinceFire +(1|Site),
              MetaRare, REML = TRUE);summary(Mod21spFinal)

Top2Mods<-tab_model(Mod17spFinal,Mod21spFinal, p.val = "kr", show.stat = TRUE);Top2Mods

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
dir.create(file.path("Analysis/MixedModels/Ectomycorrhizal"), recursive = TRUE)
capture.output(summary(MFinalSp),file="Analysis/MixedModels/Ectomycorrhizal/ASV-MixedModel-Results.csv")












