#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")

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
UnburnedSap<-MetaRareSap[which(Treatment== "Unburned"), ]
head(UnburnedSap)
dim(UnburnedSap)#102x34

# * * BURNED PLOTS ------------------------------
BurnedSap<-MetaRareSap[which(Treatment== "Burned"), ]
head(BurnedSap);dim(BurnedSap)#106x34
detach(MetaRareSap)


#
#
#**************************************************************************************************************-----
# NON METRICS MULTIDIMENSIONAL SCALING ---------------------------------------------------------------------
#**************************************************************************************************************-----
par=(mfrowc(1,1))

#Calculate distance metric ...................................................
library(vegan)
Distance2Sap<-avgdist(SaprobeTrans2, 2069, iterations=100, 
                      meanfun=median, transf= sqrt, dmethod="bray")



#CALCULATE NMDS ...............................................................-----
NMDSsap<-metaMDS(Distance2Sap, autotransform = FALSE, engine = "monoMDS", k=3, 
                 weakties =TRUE, model="global",maxit = 400, try=80,
                 trymax=100);print(NMDSsap)#Stress=0.1256035

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
                theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),                             text = element_text(size=14),        
                    axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
                    axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
                    legend.title = element_text(size = 16, colour = "black", face = "bold"),
                    legend.text = element_text(size = 15,  colour ="black"), 
                    legend.position = "right")+ 
                scale_colour_manual(values = c("#a2673f","#45877f"));NMDSSaprobe
  
#NMDS W YEAR SHAPES .........................................................
NMDSSaprobeTSF<-ggplot(data.scoresSap, aes(x = NMDS1, y = NMDS2)) + 
                  geom_point(size = 3.5, aes( colour = Treatment, shape=FireYear))+ 
                  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),                             text = element_text(size=14),        
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

#Only kept the ones that were significant after running all variables.............................
attach(MetaRareSap)
FIT2<-envfit(NMDSsap ~ Aspect + TimeSinceFire + AvgDepthOM + AvgDistanceTree + TC +
              TP + TC.TP + Slope, data=MetaRareSap, 
             permutations =9999, strata=Site);FIT2



#capture.output(FIT2, file="Analysis/Diversity/Saprobes/Tables/EnvFit-Saprobes.csv")

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
                     theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),                             text = element_text(size=14),        
                        axis.text.y = element_text(colour = "black", size = 13, face = "bold"), 
                        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
                        legend.title = element_text(size = 16, colour = "black", face = "bold"),
                       legend.text = element_text(size = 15,  colour ="black"), 
                       legend.position = "right")+ theme(plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+
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
AdonTrtTsfSap<-adonis2(SaprobeTrans2~Treatment+TimeSinceFire, 
                    permutations=9999, method="bray", strata="Plot",
                    data=MetaRareSap);AdonTrtTsfSap

capture.output(AdonTrtTsfSap,file="Analysis/Diversity/Saprobes/Tables/Adonis-TrtTSF-Saprobes.csv")



















