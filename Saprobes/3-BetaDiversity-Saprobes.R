#Spet 27, 2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")

library(vegan)

#Load Data--------------------------------------------------------------------------
MetaRareSap<- read.csv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.csv", na.strings = "N/A", header = TRUE) 
SaprobeTrans2<-read.csv("Analysis/ASV-Tables/Saprobes/SaprobeTrans2.csv", row.names=1, check.names = FALSE)
dim(SaprobeTrans2);dim(MetaRareSap)
sort(rowSums(SaprobeTrans2), decreasing = TRUE)

##
##
#*****************************************************************************************************************
# QUALITY CONTROL )----------------------------------------------------------------------------------------------
#*****************************************************************************************************************
# Convert to factor------------------------------------------------
MetaRareSap$FireYear <- as.factor(MetaRareSap$FireYear)
MetaRareSap$TimeSinceFire <- as.factor(MetaRareSap$TimeSinceFire)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRareSap)
UnburnedSap<-MetaRareSap[which(Treatment== "Unburned"), ];dim(UnburnedSap)#102x34

# * * BURNED PLOTS ------------------------------
BurnedSap<-MetaRareSap[which(Treatment== "Burned"), ];dim(BurnedSap)#106x34
detach(MetaRareSap)


#
#
#**************************************************************************************************************-----
# NON METRICS MULTIDIMENSIONAL SCALING ---------------------------------------------------------------------
#**************************************************************************************************************-----
par=(mfrowc(1,1))
library(vegan)

#Calculate distance metric ................................................................................
Distance2Sap<-avgdist(SaprobeTrans2, 2069, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")



#CALCULATE NMDS ...........................................................................................
NMDSsap<-metaMDS(Distance2Sap, autotransform = FALSE, engine = "monoMDS", k=3, weakties =TRUE,
                 model="global",maxit = 400, try=80,trymax=100);print(NMDSsap)#Stress=0.1256035

#
#
# * * GOODNESS OF FIT-TEST ---------------------------------------------------------
gof.nmds<-goodness(NMDSsap)
plot(NMDSsap, display="sites", cex=2*gof.nmds/mean(gof.nmds))#display GOF in graph
stressplot(NMDSsap, p.col= "blue", l.col="red", lwd=2)#0.984


#**********************************************************************************************************------
# NMDS PLOTS --------------------------------------------------------------------------------------------------
#*********************************************************************************************************-------

# * PLOT 1-FIREYEAR -------------------------------------------------------------------------
#extract NMDS scores (x and y coordinates)
data.scoresSap = as.data.frame(scores(NMDSsap))

#add columns to data frame, based on the values you want to plot
data.scoresSap$Treatment = MetaRareSap$Treatment
data.scoresSap$FireYear= MetaRareSap$FireYear
data.scoresSap$TimeSinceFire = MetaRareSap$TimeSinceFire


# CREATE PLOTS ...............................................................................
NMDSSaprobe<-ggplot(data.scoresSap, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3.5, aes( colour = Treatment))+ 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), text = element_text(size=14), 
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"), 
        legend.position = "right")+ 
        scale_colour_manual(values = c("#a2673f","#45877f"));NMDSSaprobe
  
#NMDS W YEAR SHAPES .........................................................
NMDSSaprobeTSF<-ggplot(data.scoresSap, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment, shape=FireYear))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"), 
        legend.position = "right")+ 
  scale_colour_manual(values = c("#a2673f","#45877f")); NMDSSaprobeTSF

#
#
#****************************************************************************************************************-----
# ----------------------------------    PLOT NMDS W ENVIFIT ARROWS    -----------------------------------------
#****************************************************************************************************************-----
#Interpretation: fit EV arrow points to the direction of most rapid change in the the EV variable.Direction is 
#proportional to the correlation between ordination and environmental variable. Often this is called the
#strength of the gradient

attach(MetaRareSap)

#Only kept the ones that were significant after running all variables...................................
FIT2<-envfit(NMDSsap ~AvgDepthOM+  TP +  TC.TP + TC + AvgDistanceTree +  Slope, data=MetaRareSap, 
             permutations =9999, strata=Site);FIT2

capture.output(FIT2, file="Analysis/Diversity/Saprobes/Tables/EnvFit-Saprobes.csv")

               
#Extract Data from envfit object for plotting ..........................................................
env.scores2<- as.data.frame(scores(FIT2, display = "vectors")) #extracts scores
env.scores2<- cbind(env.scores2, env.variables = rownames(env.scores2)) #give them names



#multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
NMDSsapTSFenv <- ggplot(data.scoresSap, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment, shape=TimeSinceFire))+ 
  geom_segment(data = env.scores2,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.20, "cm")), colour = "black") + 
  geom_text(data = env.scores2,aes(x = NMDS1, y = NMDS2, label=env.variables), 
            size = 5, hjust= 0.7, vjust= -0.5)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), text = element_text(size=14),  
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"), 
        legend.position = "right")+ 
  theme(plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+
  scale_colour_manual(values = c("#a2673f","#45877f"));NMDSsapTSFenv
detach(MetaRareSap)


#*************************************************************************************************************------
#EXPORT GRAPHS------------------------------------------------------------------------------------------------------
#*************************************************************************************************************------

pdf("Analysis/Diversity/Saprobes/Graphs/Saprobe-NMDS-Trt.pdf", height=6, width=8)
NMDSSaprobe
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/Saprobe-NMDS-Trt&TSF.pdf", height=6, width=8)
NMDSSaprobeTSF
dev.off()

pdf("Analysis/Diversity/Saprobes/Graphs/Saprobe-NMDS-Trt&TSF-w-Envfit.pdf", height=6, width=8)
NMDSsapTSFenv
dev.off()


#************************************************************************************************************I**----
#TEST FOR SIGNIFICANCE ---------------------------------------------------------------------------------------------
#***************************************************************************************************************----

#Permanova#control permutations to plot, to account for pseudoreplication
#Adonis2 using Site, Plot and or Site+Plot as random effect made no difference in overall results
# all resulted in p values of 1e-04


#*************groups have significantly different compositions (p<0.05)..................................
AdonTrtTsfSap<-adonis2(SaprobeTrans2~Treatment * TimeSinceFire + Elevation + Aspect+
                  Slope + TC.TP + TC.TN+Treatment*AvgDepthOM,strata="Site", 
                permutations=9999, method="bray", data=MetaRareSap);AdonTrtTsfSap


capture.output(AdonTrtTsfSap,file="Analysis/Diversity/Saprobes/Tables/Adonis-TrtTSF-Saprobes.csv")


#TEST OF HOMOGENEITY --------------------------------------------------------------------------------------------
#Betadisper calculates the average distance of group members to the group centroid 
#in multivariate space (generated by a distance matrix). 
#Then, an ANOVA is done to test if the dispersions (variances) of groups are different.
#("Null hypothesis of no difference in dispersion between groups";

#BETADISP...........ANOVA's p-value is not significant meaning that group dispersions are homogenous 
AVbetaTrt<-anova(betadisper(Distance2Sap, MetaRareSap$Treatment));AVbetaTrt# p=0.51 signif. different spreads---groip homogeneous  
AVbetaTSF<-anova(betadisper(Distance2Sap, MetaRareSap$TimeSinceFire));AVbetaTSF#p=5.92e-08 signif. different spreads ---group not homogeneous only 11 year in comparison to other years
AVbetaAsp<-anova(betadisper(Distance2Sap, MetaRareSap$Aspect));AVbetaAsp#7.514e-16 ***

TkbetaTrt<-TukeyHSD(betadisper(Distance2Sap, MetaRareSap$Treatment));TkbetaTrt
TkbetaTsf<-TukeyHSD(betadisper(Distance2Sap, MetaRareSap$TimeSinceFire));TkbetaTsf#only 11-2,11-3, 11-5
TkbetaAsp<-TukeyHSD(betadisper(Distance2Sap, MetaRareSap$Aspect));TkbetaAsp# S-E, SE-S, SW-SE



#Visualize distance to group centroid ..............................................................
par(mfrow=c(1,3))
boxplot(betadisper(Distance2Sap, MetaRareSap$Treatment))
boxplot(betadisper(Distance2Sap, MetaRareSap$TimeSinceFire))
boxplot(betadisper(Distance2Sap, MetaRareSap$Aspect))

#If close to one another= similar in composion
par(mfrow=c(1,1))
plot(betadisper(Distance2Sap, MetaRareSap$Treatment), hull=FALSE, ellipse=TRUE)
plot(betadisper(Distance2Sap, MetaRareSap$TimeSinceFire),hull=FALSE, ellipse=TRUE)
plot(betadisper(Distance2Sap, MetaRareSap$Aspect),hull=FALSE, ellipse=TRUE)



#Results= there is an effect of Trt on saprobic communities (adonis)
#and that the communities from each Treatment display similar homogeneity.

#Anova
capture.output(AVbetaTrt,file="Analysis/Diversity/Saprobes/Tables/Beta/Anova-beta-Trt.csv")
capture.output(AVbetaTSF,file="Analysis/Diversity/Saprobes/Tables/Beta/Anova-beta-Tsf.csv")
capture.output(AVbetaAsp,file="Analysis/Diversity/Saprobes/Tables/Beta/Anova-beta-Aspect.csv")

#Tukey
capture.output(TkbetaTrt,file="Analysis/Diversity/Saprobes/Tables/Beta/Tukey-beta-Trt.csv")
capture.output(TkbetaTsf,file="Analysis/Diversity/Saprobes/Tables/Beta/Tukey-beta-tsf.csv")
capture.output(TkbetaAsp,file="Analysis/Diversity/Saprobes/Tables/Beta/Tukey-beta-Aspect.csv")


#***************************************************************************************************-----
#------------- EXPLANATION OF RESTULTS -----------------------------------------------------------------
#***************************************************************************************************-----
#
#--------------TREATMENT .............................................................
# So our groups (Burned vs Unburned) present homogeneity among group dispersions 
# (compositions vary similarly) while having significantly different compositions.

#--------------TIME SINCE FIRE  ........................................................
# So our groups (Time since Fire: 2, 3, 5, 11) do not have homogeneity among group dispersions 
# due to comparison of group 11 to 2, 3, 5 (compositions vary), however, 
# they have  significantly different compositions.














#Full Model to test environmental variables ............................................................
AdonTrtTsfSap<-adonis2(SaprobeTrans2~Treatment*TimeSinceFire + Treatment*AvgDepthOM+ 
                TC.TP + TC.TN + TN + TC + TP+ AvgDistanceTree + Slope +
                  Elevation + Aspect, strata="Site", permutations=999, method="bray", 
                data=MetaRareSap);AdonTrtTsfSap


capture.output(AdonTrtTsfSap,file="Analysis/Diversity/Saprobes/Tables/Adonis-TrtTSF-Saprobes.csv")


















