#9-14-2020

#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")

#Load libraries-------------------------------------
library(plyr)
library(tidyverse)
library(vegan)
library(EcolUtils)#rrarefy.perm


#LOAD DATAFILES----------------------------------------------------------------------------------------------------
metadata<-read.csv("Metadata/Metadata-PIPO.csv", na.strings = "N/A", header = TRUE) #Full Metadata
RawOtu<-read.csv("Qiime/FunGuild/Fungal-table-with-taxonomy.guilds_editted.csv", row.names = 1) #Raw OTU table w/o singletons, contaminated sample removed
dim(RawOtu)#8585  223
head(RawOtu[1:3,])

#Subset data to maintain only those taxa that are classified as Saprobes 
RawSaprobes<-RawOtu[which(RawOtu$EdittedGuild =="Saprotroph"),]
dim(RawSaprobes)#2827x 224
(RawSaprobes$Guild)#All saprobes


# QUALITY CONTROL-------------------------------------------------------------------------------------------------------
# * Remove FeatureID and Taxonomy Labels to use later-------------
lastcolumn <- ncol(RawSaprobes); lastcolumn #taxonomy colum
TaxonGuild<-RawSaprobes[,214:225]; head(TaxonGuild[1:2,]) #Extract taxonomy column from the dataset
featureID <- row.names(RawSaprobes); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,TaxonGuild); TaxonID #Create a dataframe of the featureID taxonomy and guild

# * Clean Taxonomy labels by separating them by (;)--------------------------------------------------------------------------------
   #divide last row into subsets #separate taxonomy w spaces 

SplitTaxonomy<-ldply(str_split(string = TaxonID$taxonomy, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));SplitTaxonomy2
tail(SplitTaxonomy2[1:2,])
dim(SplitTaxonomy2)#2827x7

# * Change the format of the TaxonID dataframe created above so that is is split
TaxonID <- cbind(TaxonID[,1:12 ], SplitTaxonomy2);head(TaxonID[1:2,])
rownames(TaxonID) <- featureID #add feaure ID to table
tail(TaxonID)

# * * Look at your dataset to make sure that they have the same length------------------------------------------------------------
dim(TaxonID);dim(RawSaprobes) # should have same row length  (890)


# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE without ; separations
RawSaprobe2<- cbind(RawSaprobes, TaxonID)
head(RawSaprobe2)

# CREATE OTU TABLE--NO TAXONOMY --------------------------------------------------------------------------------------------------
RawSaprobeTable  <- RawSaprobe2[ ,1:213]#number as shown above in taxonguild
head(RawSaprobeTable[1:2,])
dim(RawSaprobeTable)#2827x213

# TRANSPOSE RAW OTU TABLE-----------------------------------------------------------------------------------
SaprobeTrans <- t(RawSaprobeTable) #transpose OTU table 
dim(SaprobeTrans) #213x2827

#EXPORT DATA ........................................................................
dir.create("Analysis/ASV-Tables/Saprobes")
dir.create("Analysis/QualityControl/Saprobes")

write.csv(RawSaprobeTable, "Analysis/ASV-Tables/Saprobes/RawSaprobeTable.csv")
write.csv(SaprobeTrans, "Analysis/ASV-Tables/Saprobes/RawSaprobeTrans.csv")
write.csv(TaxonID, "Analysis/QualityControl/Saprobes/TaxoSaprobes.csv")


##
##
#*******************************************************************************************************----
# QUALITY CONTROL  -----------------------------------------------------------------------------------------
#*******************************************************************************************************----
#    NegControl look looks good, no mock community was used

# * look at Negative DNA controls .............................................
NegControls <- which(colnames(RawSaprobeTable) %in% c("NC")); NegControls


#Redo RepseQ to see what you get when you run in R---they match-- do n ........
SaprobeTotSeq<-sort(rowSums(SaprobeTrans), decreasing = TRUE);SaprobeTotSeq 


##
#**********************************************************************************************************----
# RAREFY FUNGAL DATA  -----------------------------------------------------------------------------------------
#**********************************************************************************************************----

# Chose rarefaction value using qiime feature table

# * * Trial 1 Remove samples w not enoug sequences............................................
#Subset data to maintain only rows with totals above rarefaction depth
SaprobeTrans2<-SaprobeTrans[rowSums(SaprobeTrans)>2069,] 
dim(SaprobeTrans2)#208x2827

# * * Rarefy table to remove low abundance sequences..........................................
#  but keep rest, normalize to 110 seq/sample, normalize 1x----

SaprobeRare1 <- rrarefy(SaprobeTrans2,  2069) 
dim(SaprobeRare1)#208x2827

# *  * Normalize the data 100x and get the mean...........................................----
SaprobeRare2<- rrarefy.perm(SaprobeTrans2, sample =2069, n = 250, round.out = T)
dim(SaprobeRare2) # 208x2867


# * * * TEST IF RAREFACTION METHOD MATTERS ...............................................----
mantel(SaprobeRare1, SaprobeRare2)# 0.9946
MantelCorr<-plot(SaprobeRare1, SaprobeRare2)


# REMOVE COLUMNS WITH ZERO READS (OTURARE2)...............................................----
# All analysis will be performed on this table---------------
# * Create an object containing zeros -------------------
zeroes <- which(colSums(SaprobeRare2)==0)
head(zeroes)

# * Remove the columns containing zero's from OTU table .................
SaprobeRare3 <- SaprobeRare2[ , -zeroes]
head(SaprobeRare3);dim(SaprobeRare3)#208x2727

# QUICK RAREFACTION CURVES (SPECIES ACCUMULATION CURVES) ------------------------------------------------------
# Has its own code


# * * EXPORT RAREFACTION TABLES  ------------------------------------------------------------------------------
dir.create("Analysis/QualityControl/Saprobes/Graphs/")
write.csv(SaprobeRare3, "Analysis/ASV-Tables/Saprobes/RareSaprobe.csv")#table w/o zeroes
write.csv(SaprobeTrans2, "Analysis/ASV-Tables/Saprobes/RareSaprobeTrans2.csv")#table w/o zeroes
write.csv(SaprobeTotSeq, "Analysis/QualityControl/Saprobes/SaprobeRepSeq.csv")



#This code takes a long time to run.....................................................................
pdf("Analysis/QualityControl/Saprobes/MantelTest-Rarefaction.pdf", height=6, width=5)
MantelCorr
dev.off()



## Remove samples that were removed during rarefaction from metadata ....................................
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3

metadata<-read.csv("Metadata/Metadata-PIPO.csv", na.strings = "N/A", header = TRUE) #Full Metadata

# * * LOAD OTU TABLES -----------------------------------------------------------------------------------
# Reload Raw OTU table W/O Zeros (ADDED column name renamed to SampleID)
SaprobeRare3<-read.csv("Analysis/ASV-Tables/Saprobes/RareSaprobe.csv",check.names = FALSE) 
dim(SaprobeRare3)#208x2732
names(SaprobeRare3)[1]<-"SampleID"

#remove unwanted samples to rarefy matadata
MetaRareSaprobe<-metadata %>% semi_join(SaprobeRare3, by = "SampleID") # keep rows with matching ID
dim(MetaRareSaprobe)#208x26 (verity same rows as above)

#Export rarefied metadata --------------------------------------------------------------------------------
write.csv(MetaRareSaprobe, "Analysis/Metadata/Saprobes/MetaRareSaprobe.csv")

##
#**************************************************************************************************************************----
# CALCULATE ALPHADIVERITY METRICS  --------------------------------------------------------------------------------------------
#**************************************************************************************************************************----


# RICHNESS METHOD 1: BIODIVERSITY CALCULATED..........................................................................----
# * Load metadata that correspands to R-rarefied table created above -------
attach(MetaRareSaprobe)
SaprobeRare3<-read.csv("Analysis/ASV-Tables/Saprobes/RareSaprobe.csv", row.names = 1, check.names = FALSE)#reload but put SampleID as row.names
MetaRareSaprobe<-read.csv("Analysis/Metadata/Saprobes/MetaRareSaprobe.csv")

# * Richness calculations using Biodivesity package...........................................................----
SaprobeRichness <- estimateR(SaprobeRare3)
SaprobeRichness<- t(estimateR(SaprobeRare3))#run function on transformed table to switch rows and columns



# CREATE FUNCTION TO ESTIMATE DIVERSITY METRICS -------------------------------------------------------------------
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

# DIVERSITY INDICES ..................................................................................
shanEntro <- diversity(SaprobeRare3) # Shannon entropy
shannon <- exp(SaprobeRare3 ) ## Shannon number of diversity
simpson <- diversity(SaprobeRare3, "simpson")#Simpson diversity
simpEven<- diversity(SaprobeRare3, "inv") ## Simpson Inverse (Eveness)

# * Dataframe of shannon entropy, diversity &  simpson diversity--------------
Saprobe.richness <- data.frame(shanEntro, shannon, simpson, simpEven)

# * Add above diversity metrics to S obs, chao1, ACE
SaprobeSppRichness <- cbind(Saprobe.richness, SaprobeRichness)
dim(SaprobeSppRichness)#208X2739


#
#****************************************************************************************************************----
# EXPORT DIVERSITY DATA AND ATTACH TO RAREFIED METADATA ---------------------------------------------------
#****************************************************************************************************************----
dir.create("Analysis/Diversity/Saprobes")
write.csv(SaprobeSppRichness, "Analysis/Diversity/Saprobes/SaprobeCoreMetricsAll.csv") 


#sUBSET ONLY THE RICHNESS METRICS, they are at end of file and one in 1st column.................
ncol(SaprobeSppRichness)# 2739

#1 column and 6th column
CoreMetrics1<-SaprobeSppRichness[,2733:2739];head(CoreMetrics1)#2739-6
CoreMetrics2<-SaprobeSppRichness[,1];head(CoreMetrics2)

CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)
names(CoreMetrics)[7]<- "ShannonEntropy"


#Export core metrics results----------------------------------------------
write.csv(CoreMetrics, "Analysis/Diversity/Saprobes/CoreMetrics.csv")

#****************************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare ---------------------------------------------------
#****************************************************************************************************************----

CoreMetrics<-read.csv("Analysis/Diversity/Saprobes/CoreMetrics.csv", check.names = FALSE)
names(CoreMetrics)[1]<- "SampleID"


MetaRareSpp<-merge(MetaRareSaprobe,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#208X33

#Export rarefied metadata w core metrics attached---------------------
dir.create(file.path("Analysis/Metadata/Saprobes/"))

write.csv(MetaRareSpp, "Analysis/Metadata/Saprobes/MetaRareCoreSaprobeNew.csv")

