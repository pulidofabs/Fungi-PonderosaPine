#9-15-2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")

#LOAD PACKAGES ------------------------------------------------------------------------------------------
library(tidyverse) #required to load tax table

#Plotting  
library(ggplot2)
library(ggpubr) 
library(tidyr)#manipulate data
library(ggfortify) #Manipulate and graph plots
library(ochRe)# color palette library(multcomp)#multiple comparisons--for glm 


#Analysis
library(nlme) #nonlinear models
library(multcomp)
library("stringr")
library("plyr")

library(EcolUtils) #install.packages("devtools");library(devtools); devtools::install_github("GuillemSalazar/EcolUtils")
library(SPECIES)
library(BiodiversityR)
library(scales)

#LOAD DATAFILES----------------------------------------------------------------------------------------------------
# * * LOAD METADATA --------------------------------------------------------
metadata<-read.csv("Metadata/Metadata-PIPO.csv", na.strings = "N/A", header = TRUE) #Full Metadata

# * * LOAD OTU TABLES -----------------------------------------------------
RawOtu<-read.csv("Qiime/FunGuild/Fungal-table-with-taxonomy.guilds_editted.csv", row.names = 1) #Raw OTU table w/o singletons, contaminated sample removed
dim(RawOtu)#8585  223


# QUALITY CONTROL-------------------------------------------------------------------------------------------------------
# * Remove FeatureID and Taxonomy Labels to use later-------------
lastcolumn <- ncol(RawOtu); lastcolumn #taxonomy colum
TaxonGuild<-RawOtu[,214-223]; head(TaxonGuild) #Extract taxonomy column from the dataset
featureID <- row.names(RawOtu); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,TaxonGuild); TaxonID #Create a dataframe of the featureID and taxonomy

# * Clean Taxonomy labels by separating them by (;)--------------------------------------------------------------------------------
   #divide last row into subsets #separate taxonomy w spaces 

TaxonID$taxonomy
SplitTaxonomy<-ldply(str_split(string = TaxonID$taxonomy, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));SplitTaxonomy2
tail(SplitTaxonomy2)
dim(SplitTaxonomy2)

# * Change the format of the TaxonID dataframe created above so that is is split
TaxonID <- cbind(TaxonID[,1:7 ], SplitTaxonomy2 );TaxonID
rownames(TaxonID) <- featureID #add feaure ID to table
tail(TaxonID)

# * * Look at your dataset to make sure that they have the same length------------------------------------------------------------
dim(TaxonID);dim(RawOtu) # should have same row length 


# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE
RawOtu2<- cbind(RawOtu, TaxonID)
head(RawOtu2)

# CREATE OTU TABLE--NO TAXONOMY --------------------------------------------------------------------------------------------------
RawOtuTable  <- RawOtu2[ ,1:213]#number as shown above-1
head(RawOtuTable)

# TRANSPOSE RAW OTU TABLE---------------------------------------------------------------------------------------------------------
OtuTrans <- t(RawOtuTable) #transpose OTU table 
Dim(OtuTrans) #213-8585

#EXPORT DATA ----------------------------------------dir.create(file.path("Fungi/Output/Abundance/Tables"))
dir.create("Analysis/ASV-Tables")
dir.create("Analysis/QualityControl")

write.csv(RawOtuTable, "Analysis/ASV-Tables/RawOtuTable.csv")
write.csv(OtuTrans, "Analysis/ASV-Tables/OtuTrans.csv")
write.csv(TaxonID, "Analysis/QualityControl/TaxononomyLabels.csv")


#
#
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# QUALITY CONTROL  ----------------------------------------------------------------------------
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*
#    NegControl look looks good, no mock community was used

# * look at Negative DNA controls ----------------------
NegControls <- which(colnames(RawOtuTable) %in% c("NC")); NegControls


#Redo RepseQ to see what you get when you run in R---they match-- do n
TotSeq<-sort(rowSums(OtuTrans), decreasing = TRUE);TotSeq#will not export it is the same

#.----
#.----
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*------------
# RAREFY TABLE- FULL FUNGAL--------------------------------------------------------------------------------------------------------
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*------------
# Chose rarefaction value using qiime feature table

# * * Trial 1 Remove samples w not enoug sequences----------------------------------------------------------------------
OtuTrans2<-OtuTrans[rowSums(OtuTrans)>6229,] #Subset data to maintain only rows with totals above rarefaction depth
dim(OtuTrans2)#210X8585

# * * Rarefy table to remove low abundance sequences but keep rest, normalize to 110 seq/sample, normalize 1x----
OtuRare1 <- rrarefy(OtuTrans2,  6229) #not sure what is the difference btwn this one and other one
dim(OtuRare1)#210X8585

# *  * Normalize the data 100x and get the mean------------------------------------
OtuRare2<- rrarefy.perm(OtuTrans2, sample =6229, n = 250, round.out = T)
dim(OtuRare2) # 210x8585


# * * * TEST IF RAREFACTION METHOD MATTERS ----------------------------------------
mantel(OtuRare1, OtuRare2)# 0.9997
MantelCorr<-plot(OtuRare1, OtuRare2)


# REMOVE COLUMNS WITH ZERO READS (OTURARE2)----------------------------------------
# * Create an object containing zeros -------------------
zeroes <- which(colSums(OtuRare2)==0)
head(zeroes)

# * Remove the columns containing zero's from OTU table ----------------------------
OtuRare3 <- OtuRare2[ , -zeroes]
head(OtuRare3)
dim(OtuRare3)#8464

# RAREFACTION CURVES (SPECIES ACCUMULATION CURVES) ----------------------------------



# * * EXPORT RAREFACTION TABLES  -----------------------------------------------------------------
dir.create("Analysis/QualityControl")
dir.create("Analysis/QualityControl/Graphs")
write.csv(OtuTrans2, "Analysis/ASV-Tables/OtuTrans2.csv")#
write.csv(OtuRare2, "Analysis/ASV-Tables/OtuRare.csv")#table w/o zeroes
write.csv(OtuRare3, "Analysis/ASV-Tables/OtuRare3-NZ.csv")#table w/o zeroes

pdf("Analysis/QualityControl/Graphs/MantelTest-Rarefaction.pdf", height=6, width=5)
MantelCorr
dev.off()



## Remove samples that were removed during rarefaction from metadata 
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3
metadata<-read.csv("Metadata/Metadata-PIPO.csv", na.strings = "N/A", header = TRUE) #Full Metadata

# * * LOAD OTU TABLES -----------------------------------------------------
OtuRare3<-read.csv("Analysis/ASV-Tables/OtuRare.csv") #Raw OTU table w/o singletons, contaminated sample removed
dim(OtuRare3)#210x8585

#remove unwanted samples to rarefy matadata
MetaRare<-metadata %>% semi_join(OtuRare3, by = "SampleID") # keep rows with matching ID


#Export rarefied metadata ------------------
write.csv(MetaRare, "Analysis/Metadata/Fungi/MetaRare.csv")




#.----
#.----
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--------------
# CALCULATE ALPHA DIVERSITY --------------------------------------------------------------------------------------------------------
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--------------

# RICHNESS METHOD 1: BIODIVERSITY CALCULATED------------------------------------------------------------------------------
# * Load metadata that correspands to R-rarefied table created above -------

OtuRare3<-read.csv("Analysis/ASV-Tables/OtuRare.csv", row.names = 1)
MetaRare<-read.csv("Analysis/Metadata/Fungi/MetaRare.csv")
attach(MetaRare)

# * Richness calculations using Biodivesity package-----------------------------------
OtuRichness <- estimateR(OtuRare3)
OtuRichness<- t(estimateR(OtuRare3))#run function on transformed table to switch rows and columns

# CREATE FUNCTION TO ESTIMATE DIVERSITY METRICS ------------------------------
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  # pdf(paste("figures/richnesscores_",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  mtext("B",side=3,adj=0)
  #dev.off()
}

# DIVERSITY INDICES -------------------------------------------------------------------------------
shanEntro <- diversity(OtuRare3) # Shannon entropy
shannon <- exp(OtuRare3 ) ## Shannon number of diversity
simpson <- diversity(OtuRare3, "simpson")#Simpson diversity
simpEven<- diversity(OtuRare3, "inv") ## Simpson Inverse (Eveness)

# * Dataframe of shannon entropy, diversity &  simpson diversity--------------
otu.richness <- data.frame(shanEntro, shannon, simpson, simpEven)

# * Add above diversity metrics to S obs, chao1, ACE
OtuSppRichness <- cbind(otu.richness, OtuRichness)
dim(OtuSppRichness)#210x8593


#****************************************************************************************************************----
# EXPORT DIVERSITY DATA AND ATTACH TO RAREFIED METADATA ---------------------------------------------------
#****************************************************************************************************************----
# * * * Export Biodiversity package spp richness--------------------------------------------------------------------
dir.create("Analysis/Diversity/Fungi")
write.csv(OtuSppRichness, "Analysis/Diversity/Fungi/EcMCoreMetricsAll.csv") #all above species calculations

#sUBSET ONLY THE RICHNESS METRICS, they are at end of file and one in 1st column.................
ncol(OtuSppRichness)#8593

#1 column and 6th column
CoreMetrics1<-OtuSppRichness[,8587:8593];head(CoreMetrics1)#8593-6 
CoreMetrics2<-OtuSppRichness[,1];head(CoreMetrics2)

CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)
names(CoreMetrics)[7]<- "ShannonEntropy"; head(CoreMetrics)


#Export core metrics results----------------------------------------------
dir.create("Analysis/Diversity/Fungi")
write.csv(CoreMetrics, "Analysis/Diversity/Fungi/CoreMetrics.csv")

##
#****************************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare ---------------------------------------------------
#****************************************************************************************************************----

CoreMetrics<-read.csv("Analysis/Diversity/Fungi/CoreMetrics.csv", check.names = FALSE)
names(CoreMetrics)[1]<- "SampleID"


MetaRareSpp<-merge(MetaRare,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#210x33

#Export rarefied metadata w core metrics attached------------------------------------------
write.csv(MetaRareSpp, "Analysis/Metadata/Fungi/MetaRareCoreFun.csv")




