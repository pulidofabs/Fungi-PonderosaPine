#Last date updated
#9/27

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")



library(ggplot2)
library(vegan)


#Load Data--------------------------------------------------------------------------
MetaRareEcM<- read.csv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv", na.strings = "N/A", header = TRUE) 
EcMTrans2<-read.csv("Analysis/ASV-Tables/Ectomycorrhizal/EcMTrans2.csv", row.names=1, check.names = FALSE)



#*****************************************************************************************************************
# QUALITY CONTROL )----------------------------------------------------------------------------------------------
#*****************************************************************************************************************

# Convert to factor------------------------------------------------
MetaRareEcM$FireYear <- as.factor(MetaRareEcM$FireYear)
MetaRareEcM$TimeSinceFire <- as.factor(MetaRareEcM$TimeSinceFire)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRareEcM)
Unburned<-MetaRareEcM[which(Treatment== "Unburned"), ]
head(Unburned)
dim(Unburned)#105x33

# * * BURNED PLOTS ------------------------------
Burned<-MetaRareEcM[which(Treatment== "Burned"), ]
head(Burned)
dim(Burned)#87x33
detach(MetaRareEcM)

#
#
#*******************************************************************************************************************-----
# NON METRICS MULTIDIMENSIONAL SCALING ----------------------------------------------------------------------------------
#*******************************************************************************************************************_----
library(vegan)

# CALCULATE DISTANCE.......................................................
Distance2<-avgdist(EcMTrans2, 118, iterations=100, meanfun=median,
                   transf= sqrt, dmethod="bray" )

#CALCULATE NMDS ............................................................
NMDS<-metaMDS(Distance2, autotransform = FALSE, engine = "monoMDS", k=3, 
              weakties =TRUE, model="global",maxit = 400, try=80,
              trymax=100);print(NMDS)#Stress=0.1482633



# * * GOODNESS OF FIT TEST ..................................................
par(mfrow=c(1,1))
gof.nmds<-goodness(NMDS)
GOFplot<-plot(NMDS, display="sites", cex=2*gof.nmds/mean(gof.nmds))
StressPlot<-stressplot(NMDS, p.col= "blue", l.col="red", lwd=2)#

#
#**********************************************************************************************----
#-------------------------------- NMDS PLOTS ------------------------------------------------
#**********************************************************************************************----
#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(NMDS))

#add columns to data frame, based on the values you want to plot
data.scores$Treatment = MetaRareEcM$Treatment
data.scores$FireYear= MetaRareEcM$FireYear
data.scores$TimeSinceFire = MetaRareEcM$TimeSinceFire
data.scores$Aspect = MetaRareEcM$Aspect#no difference, all clumped
data.scores$Elevation = MetaRareEcM$Elevation# 1100-1162 maybe
data.scores$AvgDistanceTree = MetaRareEcM$AvgDistanceTree# NO
data.scores$TP = MetaRareEcM$TP#YES
data.scores$AvgDepthOM = MetaRareEcM$AvgDepthOM#YES


#PLOT ..................................................................................
NMDSecm<-ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"), 
        legend.position = "right")+
        scale_colour_manual(values = c("#a2673f","#45877f"));NMDSecm


NMDSecmTSF<-ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3.5, aes( colour = Treatment, shape=TimeSinceFire))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),                             
        text = element_text(size=14),        
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"),
        legend.position = "right")+ 
        scale_colour_manual(values = c("#a2673f","#45877f"));NMDSecmTSF

#
#
#****************************************************************************************************************-----
# ----------------------------------    PLOT NMDS W ENVIFIT ARROWS    -----------------------------------------
#****************************************************************************************************************-----
#Interpretation: fit EV arrow points to the direction of most rapid change in the the EV variable.Direction is 
#proportional to the correlation between ordination and environmental variable. Often this is called the
#strength of the gradient
#Only kept the ones that were significant after running all variables.........................

attach(MetaRareEcM)
FIT<-envfit(NMDS ~ AvgDepthOM + Elevation, data=MetaRareEcM, permutations =9999, na.rm = TRUE, strata=Site);FIT



#Export envfit results........................................................................
capture.output(FIT, file="Analysis/Diversity/Ectomycorrhizal/Tables/Envfit-EcM.csv")


# data for the envfit arrows..................................................................
env.scores<- as.data.frame(scores(FIT, display = "vectors")) #extracts scores 
env.scores<- cbind(env.scores, env.variables = rownames(env.scores))#and names to scores



#multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
NMDSecmTSFenv <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 3.5, aes( colour = Treatment, shape=TimeSinceFire))+ 
        geom_segment(data = env.scores,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                     arrow = arrow(length = unit(0.20, "cm")), colour = "black") + 
        geom_text(data = env.scores,aes(x = NMDS1, y = NMDS2, 
                label=env.variables),size = 5, hjust= 0.7, vjust= -0.5)+
        theme_bw()+ 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              text = element_text(size=14),
              axis.text.y = element_text(colour = "black",  size = 13, face = "bold"),
              axis.text.x = element_text(colour = "black",size = 13, face = "bold"), 
              legend.title = element_text(size = 16, colour = "black", face = "bold"),
              legend.text = element_text(size = 15,  colour ="black"), 
              legend.position = "right")+ theme(plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+
              scale_colour_manual(values = c("#a2673f","#45877f"));NMDSecmTSFenv
detach(MetaRareEcM)



#*************************************************************************************************************------
#EXPORT GRAPHS------------------------------------------------------------------------------------------------------
#*************************************************************************************************************------

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/EcM-NMDS-Trt.pdf", height=6, width=8)
NMDSecm
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/EcM-NMDS-Trt-TSF.pdf", height=6, width=8)
NMDSecmTSF
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/EcM-NMDS-Envfit.pdf", height=6, width=8)
NMDSecmTSFenv
dev.off()


#
#
#*******************************************************************************************************************-----
# --------------------STATISTICAL ANALYSIS  -----------------------------------------------------------------------------
#*******************************************************************************************************************_----
#Nested Adonis in what package .......................................................................
attach(MetaRareEcM)
library(primer)


#Permanova#control permutations to plot or site, to account for pseudoreplication
### Example of use with strata, for nested (e.g., block) designs.

#Ran a full permanova to determine which variables explain most of the variance
AdonTrtTsfEcM1<-adonis2(EcMTrans2 ~ Treatment * TimeSinceFire + Aspect + 
                  Elevation + Slope + TC.TP +  Treatment*AvgDepthOM,
                  strata=Site, perm=9999,method="bray", 
                  data=MetaRareEcM);AdonTrtTsfEcM1

capture.output(AdonTrtTsfEcM1,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/1-Adonis-All-Ectomycorrhizal.csv")


#TEST OF HOMOGENEITY --------------------------------------------------------------------------------------------
#Betadisper calculates the average distance of group members to the group centroid 
#in multivariate space (generated by a distance matrix). 
#Then, an ANOVA is done to test if the dispersions (variances) of groups are different.
#("Null hypothesis of no difference in dispersion between groups";

#BETADISP...........ANOVA's p-value is not significant meaning that group dispersions are homogenous 
AVbetaTrt<-anova(betadisper(Distance2, MetaRareEcM$Treatment)); AVbetaTrt# p=0.009695 ** signif. different spreads---group not homogeneous  
AVbetaTSF<-anova(betadisper(Distance2, MetaRareEcM$TimeSinceFire));AVbetaTSF#p=0.00365 ** * signif. different spreads---group not homogeneous  
AVbetaAsp<-anova(betadisper(Distance2, MetaRareEcM$Aspect));AVbetaAsp#6.564e-10 *** signif. different spreads---group not homogeneous  
AVbetaElev<-anova(betadisper(Distance2, MetaRareEcM$Elevation));AVbetaElev#5.133e-09 ***


TkbetaTrt<-TukeyHSD(betadisper(Distance2, MetaRareEcM$Treatment));TkbetaTrt
TkbetaTsf<-TukeyHSD(betadisper(Distance2, MetaRareEcM$TimeSinceFire));TkbetaTsf#only 11-3, 11-5
TkbetaAsp<-TukeyHSD(betadisper(Distance2, MetaRareEcM$Aspect));TkbetaAsp# S-E, SE-S, SW-S


#Visualize distance to group centroid ..............................................................
par(mfrow=c(1,3))
boxplot(betadisper(Distance2, MetaRareEcM$Treatment))
boxplot(betadisper(Distance2, MetaRareEcM$TimeSinceFire))
boxplot(betadisper(Distance2, MetaRareEcM$Aspect))

#If close to one another= similar in composion
par(mfrow=c(1,1))
plot(betadisper(Distance2, MetaRareEcM$Treatment), hull=FALSE, ellipse=TRUE)
plot(betadisper(Distance2, MetaRareEcM$TimeSinceFire),hull=FALSE, ellipse=TRUE)
plot(betadisper(Distance2, MetaRareEcM$Aspect),hull=FALSE, ellipse=TRUE)


#CONCLUSION.........................................................
#Results= there is an effect of Trt on EcM communities (adonis)
#and that the communities from each Treatment display 
#different homogeneity (sign Pvalue)

#Anova
capture.output(AVbetaTrt,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Anova-beta-Trt.csv")
capture.output(AVbetaTSF,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Anova-beta-Tsf.csv")
capture.output(AVbetaAsp,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Anova-beta-Aspect.csv")

#Tukey
capture.output(TkbetaTrt,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Tukey-beta-Trt.csv")
capture.output(TkbetaTsf,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Tukey-beta-tsf.csv")
capture.output(TkbetaAsp,file="Analysis/Diversity/Ectomycorrhizal/Tables/Beta/Tukey-beta-Aspect.csv")




#***************************************************************************************************-----
#------------- EXPLANATION OF RESTULTS -----------------------------------------------------------------
#***************************************************************************************************-----
#***** ...... significant P-value 


#--------------TREATMENT ..........................................................................
# So our groups (Burned vs Unburned) present different homogeneity among group dispersions 
# (compositions vary similarly) while having significantly different compositions.

#--------------TIME SINCE FIRE  ........................................................
# So our groups (Time since Fire: 2, 3, 5, 11) different  homogeneity among group dispersions 
# due to comparison of group 11/ and 11/5 (compositions vary), however, 
# they have  significantly different compositions.



























#Full permanova ------------------NEED TO CHECK INTO THIS TO SEE IF I NEED TO DO THIS ???????????????????????????????/
AdonAll<-adonis2(EcMTrans2 ~ Treatment * TimeSinceFire+ Treatment*AvgDepthOM + Treatment*TC.TN + Treatment*TC.TP + 
                   Treatment*TC + Treatment*TP + Treatment*TN + Slope + Elevation + Aspect,
                   perm=999,method="bray", data=MetaRareEcM);AdonAll




#Tested other nested variables to see if It made a difference, but it did not 
#AdonTrtTsfEcM2<-adonis2(EcMTrans2 ~ Treatment * TimeSinceFire, strata="Site", perm=9999,method="bray", data=MetaRareEcM) ;AdonTrtTsfEcM2
#AdonTrtTsfEcM3<-adonis2(EcMTrans2 ~ Treatment * TimeSinceFire, strata="Plot", perm=9999,method="bray", data=MetaRareEcM) ;AdonTrtTsfEcM3
#AdonTrtTsfEcM4<-adonis2(EcMTrans2 ~ Treatment * TimeSinceFire, strata="Subplot", perm=9999,method="bray", data=MetaRareEcM) ;AdonTrtTsfEcM4
#AdonTrtTsfEcM<-adonis2(EcMTrans2~Treatment*TimeSinceFire, permutations=9999,  method="bray", strata="Plot", data=MetaRareEcM) ;AdonTrtTsfEcM


#Export results........................................................................................
capture.output(AdonTrtTsfEcM1,file="Analysis/Diversity/Ectomycorrhizal/Tables/Adonis-TrtTSF-SiteRandom.csv")









