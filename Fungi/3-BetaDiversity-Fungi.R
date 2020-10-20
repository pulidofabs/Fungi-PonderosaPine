#9-15-2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")

#Load libraries .......................................
library(ggplot2)
library(vegan)


#Load Data-------------------------------------------------------------------------------------------
MetaRareFun<- read.csv("Analysis/Metadata/Fungi/MetaRareCore.csv", na.strings = "N/A", header = TRUE) 
OtuTrans2<-read.csv("Analysis/ASV-Tables/OtuTrans2.csv", row.names=1)

dim(MetaRareFun)
str(OtuTrans2)


#*******************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#*******************************************************************************************************----
# Convert to factor------------------------------------------------
MetaRareFun$FireYear <- as.factor(MetaRareFun$FireYear)
MetaRareFun$TimeSinceFire <- as.factor(MetaRareFun$TimeSinceFire)


# SUBSET DATA -----------------------------------------------------
# * * UNBURNED PLOTS ------------------------------
attach(MetaRareFun)
Unburned<-MetaRareFun[which(Treatment=="Unburned"), ]
head(Unburned)
dim(Unburned)#106X33

# * * BURNED PLOTS ------------------------------
Burned<-MetaRareFun[which(Treatment== "Burned"), ]
head(Burned)
dim(Burned)#104x33
detach(MetaRareFun)


#
#
#*******************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#*******************************************************************************************************----

Distance2<-avgdist(OtuTrans2, 6229, iterations=100, meanfun=median, transf= sqrt, dmethod="bray" )


#CALCULATE NMDS........................................................................................
NMDS<-metaMDS(Distance2, autotransform = FALSE, engine = "monoMDS", k=3, weakties =TRUE, 
              model="global", maxit = 400, try=80, trymax=100);print(NMDS)#Stress=0.1425


#
# * * Goodnes of fit-test -----------------------------
par(mfrow=c(1,1))
gof.nmds<-goodness(NMDS)
plot(NMDS, display="sites", cex=2*gof.nmds/mean(gof.nmds))#display GOF in graph
stressplot(NMDS, p.col= "blue", l.col="red", lwd=2)#0.979, 0.862



#*******************************************************************************************************----
# NON METRIC MULTIDIMENSIONAL SCALING (NMDS) ---------------------------------------------------------------
#*******************************************************************************************************----
# * PLOT 1-FIREYEAR .....................................................................
MetaRareFun$TimeSinceFire<-as.factor(MetaRareFun$TimeSinceFire)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(NMDS))

#add columns to data frame, based on the values you want to plot
data.scores$Treatment = MetaRareFun$Treatment
data.scores$FireYear= MetaRareFun$FireYear
data.scores$TimeSinceFire = MetaRareFun$TimeSinceFire



#PLOT NMDS---------------------------------------------------------------------------------------------------------

NMDSfungi<-ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment))+ 
  theme(axis.text.y = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")  + 
  scale_colour_manual(values = c("#a2673f","#45877f")); NMDSfungi



NMDSfungiTSF<-ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment, shape=TimeSinceFire))+ 
  theme(axis.text.y = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")  + 
  scale_colour_manual(values = c("#a2673f","#45877f")); NMDSfungiTSF


#
#
#****************************************************************************************************************-----
# ----------------------------------    PLOT NMDS W ENVIFIT ARROWS    -----------------------------------------
#****************************************************************************************************************-----
#Interpretation: fit EV arrow points to the direction of most rapid change in the the EV variable.Direction is 
#proportional to the correlation between ordination and environmental variable. Often this is called the
#strength of the gradient
#Only kept the ones that were significant after running all variables.........................

attach(MetaRareFun)
FIT<-envfit(NMDS ~  AvgDepthOM +  TC.TP + TP + TC + AvgDistanceTree +Slope, strata = Site,
            data=MetaRareFun, permutations =9999);FIT

capture.output(FIT, file="Analysis/Diversity/Fungi/Tables/Envfit-Fungi.csv")


# data for the envfit arrows..................................................................
env.scores<- as.data.frame(scores(FIT, display = "vectors")) #extracts scores 
env.scores<- cbind(env.scores, env.variables = rownames(env.scores))#and names to scores



#multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
NMDSenvfitFun <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes( colour = Treatment, shape=TimeSinceFire))+ 
  geom_segment(data = env.scores,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.20, "cm")), colour = "black") + 
  geom_text(data = env.scores,aes(x = NMDS1, y = NMDS2, label=env.variables),
            size = 5, hjust= 0.7, vjust= -0.5)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size=14),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 15,  colour ="black"), 
        legend.position = "right")+ 
  theme(plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+
  scale_colour_manual(values = c("#a2673f","#45877f"));NMDSenvfitFun
detach(MetaRareFunEcM)




#EXPORT ALL GRAPHS ------------------------------------------------------------------------ 
pdf("Analysis/Diversity/Fungi/Graphs/NMDS-Fungi-Trt.pdf", height=6, width=8)
NMDSfungi
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/NMDS-Fungi-Trt-TSF.pdf", height=6, width=8)
NMDSfungiTSF
dev.off()


pdf("Analysis/Diversity/Fungi/Graphs/NMDS-Fungi-Envfit-Trt-TSF.pdf", height=6, width=8)
NMDSenvfitFun
dev.off()


#
#
#*******************************************************************************************************************-----
# --------------------STATISTICAL ANALYSIS  -----------------------------------------------------------------------------
#*******************************************************************************************************************_----
#Nested Adonis in what package .......................................................................


#Permanova#control permutations to plot or site, to account for pseudoreplication
### Example of use with strata, for nested (e.g., block) designs.
attach(MetaRareFun)


#Tested both Site and Plot and they both resulted in same P-values for Trt and TSF
AdonTrtTsf<-adonis2(OtuTrans2 ~ Treatment * TimeSinceFire, strata=Site, 
                        perm=9999,method="bray", data=MetaRareFun);AdonTrtTsf

capture.output(AdonTrtTsf,file="Analysis/Diversity/Saprobes/Tables/Adonis-TrtTSF-Fun.csv")




#TEST OF HOMOGENEITY --------------------------------------------------------------------------------------------
#Betadisper calculates the average distance of group members to the group centroid 
#in multivariate space (generated by a distance matrix). 
#Then, an ANOVA is done to test if the dispersions (variances) of groups are different.
#("Null hypothesis of no difference in dispersion between groups";

#BETADISP...........ANOVA's p-value is not significant meaning that group dispersions are homogenous 
AVbetaTrt<-anova(betadisper(Distance2, MetaRareFun$Treatment)); AVbetaTrt# p=9.9e-07** signif. different spreads---group not homogeneous  
AVbetaTSF<-anova(betadisper(Distance2, MetaRareFun$TimeSinceFire)); AVbetaTSF#p=1.2e-10 ** signif. different spreads---group not homogeneous  

TkbetaTrt<-TukeyHSD(betadisper(Distance2, MetaRareFun$Treatment));TkbetaTrt
TkbetaTsf<-TukeyHSD(betadisper(Distance2, MetaRareFun$TimeSinceFire));TkbetaTsf
                        #All significant, except 5-3


#Results= there is an effect of Trt on Fun communities (adonis)
#and that the communities from each Treatment display different homogeneity.


capture.output(AVbetaTrt,file="Analysis/Diversity/Fungi/Tables/Anova-beta-Trt.csv")
capture.output(AVbetaTSF,file="Analysis/Diversity/Fungi/Tables/Anova-beta-Tsf.csv")
capture.output(TkbetaTrt,file="Analysis/Diversity/Fungi/Tables/Tukey-beta-Trt.csv")
capture.output(TkbetaTsf,file="Analysis/Diversity/Fungi/Tables/Tukey-beta-tsf.csv")




#***************************************************************************************************-----
#------------- EXPLANATION OF RESTULTS -----------------------------------------------------------------
#***************************************************************************************************-----
#
#--------------TREATMENT ..........................................................................
# So our groups (Burned vs Unburned) present different homogeneity among group dispersions 
# (compositions vary similarly) while having significantly different compositions.

#--------------TIME SINCE FIRE  ........................................................
# So our groups (Time since Fire: 2, 3, 5, 11) do not have homogeneity among group dispersions 
# due to comparison of all groups except 5-3(compositions vary), however, 
# they have  significantly different compositions.













