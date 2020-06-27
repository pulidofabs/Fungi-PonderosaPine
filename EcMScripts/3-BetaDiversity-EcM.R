#Last date updated
#4/12/2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")


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
FIT<-envfit(NMDS ~ AvgDepthOM +  AvgDistanceTree + TP + TC.TP, strata = Site, 
            data=MetaRareEcM, permutations =9999);FIT

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
        geom_text(data = env.scores,aes(x = NMDS1, y = NMDS2, label=env.variables),
                  size = 5, hjust= 0.7, vjust= -0.5)+
        theme_bw()+ 
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
              text = element_text(size=14),
              axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
              axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
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
# STATISTICAL ANALYSIS  ----------------------------------------------------------------------------------
#*******************************************************************************************************************_----

#Permanova#control permutations to plot, to account for pseudoreplication

AdonTrtTsfEcM<-adonis2(EcMTrans2~Treatment+TimeSinceFire, permutations=9999, 
                    method="bray", strata="Plot", data=MetaRareEcM) ;AdonTrtTsfEcM

#Export results........................................................................................
capture.output(AdonTrtTsf,file="Analysis/Diversity/Ectomycorrhizal/Tables/Adonis-TrtTSF.csv")









