#Reset R's Brain
rm(list=ls())

# SET WORKING DIRECTORY----------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/HolyFire/HolyFire-Year1/Fungi/")

#Load libraries-------------------------------------
library(plyr)
library(tidyverse)
library(vegan)
library(EcolUtils)#rrarefy.perm



#LOAD DATAFILES-----------------------------------------------------------------------------------------
metadata<-read.csv("Analysis/Metadata/MetaRare.csv",, na.strings = "N/A", header = TRUE, row.names = 1 )
EcMRare<-read.csv("Analysis/ASV-Tables/EcM/EcM-rare.csv", header = TRUE, row.names = 1)






#******************************************************************************************************************----
#----QUALITY CONTROL--SELECTION OF FUNGAL GUILDS ----------------------------------------------------------------------
#******************************************************************************************************************----

# * Remove FeatureID and Taxonomy Labels to use later-------------
lastcolumn <- ncol(RawEcM); lastcolumn #taxonomy colum
TaxonGuild<-RawEcM[,214:223]; head(TaxonGuild[1:3,],) #Extract taxonomy column from the dataset
featureID <- row.names(RawEcM); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,TaxonGuild); TaxonID #Create a dataframe of the featureID taxonomy and guild

# * Clean Taxonomy labels by separating them by (;)--------------------------------------------------------------------------------
   #divide last row into subsets #separate taxonomy w spaces 

SplitTaxonomy<-ldply(str_split(string = TaxonID$taxonomy, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));head(SplitTaxonomy2[1:2,])
dim(SplitTaxonomy2)#889x7

# * Change the format of the TaxonID dataframe created above so that is is split
TaxonID <- cbind(TaxonID[,1:11 ], SplitTaxonomy2);head(TaxonID[1:2,])
rownames(TaxonID) <- featureID #add feaure ID to table
tail(TaxonID)

# * * Look at your dataset to make sure that they have the same length------------------------------------------------------------
dim(TaxonID);dim(RawEcM) # should have same row length  (889)


# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE without ; separations
RawEcM2<- cbind(RawEcM, TaxonID)
head(RawEcM2[1:2,])

# CREATE OTU TABLE--NO TAXONOMY -------------------------------------------------------------------------------
RawEcMTable  <- RawEcM2[ ,1:213]#number as shown above in taxonguild -1
head(RawEcMTable[1:2,])
dim(RawEcMTable)#889x213

# TRANSPOSE RAW OTU TABLE--------------------------------------------------------------------------------------
EcMTrans <- t(RawEcMTable) #transpose OTU table 
dim(EcMTrans) #213x889

#EXPORT DATA ----------------------------------------dir.create(file.path("Fungi/Output/Abundance/Tables"))
dir.create("Analysis/ASV-Tables/Ectomycorrhizal")
dir.create("Analysis/QualityControl/Ectomycorrhizal")

write.csv(RawEcMTable, "Analysis/ASV-Tables/Ectomycorrhizal/RawEcMTable.csv")
write.csv(EcMTrans, "Analysis/ASV-Tables/Ectomycorrhizal/RawEcMTrans.csv")
write.csv(TaxonID, "Analysis/QualityControl/Ectomycorrhizal/TaxononomyEcM.csv")


#
#
#**********************************************************************************************************----
# QUALITY CONTROL  --------------------------------------------------------------------------------------------
#**********************************************************************************************************----
#    NegControl look looks good, no mock community was used

# * look at Negative DNA controls ----------------------
NegControls <- which(colnames(RawEcMTable) %in% c("NC")); NegControls


#Redo RepseQ to see what you get when you run in R---they match-- do n
EcMTotSeq<-sort(rowSums(EcMTrans), decreasing = TRUE);EcMTotSeq #rarefaction depth will be really small


#
#**********************************************************************************************************----
# RAREFY FUNGAL DATA  -----------------------------------------------------------------------------------------
#**********************************************************************************************************----

# Chose rarefaction value using qiime feature table

# * * Trial 1 Remove samples w not enoug sequences----------------------------------------------------------------------
EcMTrans2<-EcMTrans[rowSums(EcMTrans)>118,] #Subset data to maintain only rows with totals above rarefaction depth
dim(EcMTrans2)#192X889

# * * Rarefy table to remove low abundance sequences but keep rest, normalize to 110 seq/sample, normalize 1x----
EcMRare1 <- rrarefy(EcMTrans2,  118) #not sure what is the difference btwn this one and other one
dim(EcMRare1)#192x890

# *  * Normalize the data 100x and get the mean------------------------------------
# #Willl use this one to move forward
EcMRare2<- rrarefy.perm(EcMTrans2, sample =118, n = 250, round.out = T)
dim(EcMRare2) # 192X890



# * * * TEST IF RAREFACTION METHOD MATTERS ----------
mantel(EcMRare1, EcMRare2)# 0.9971
MantelCorr<-plot(EcMRare1, EcMRare2);MantelCorr


# REMOVE COLUMNS WITH ZERO READS (OTURARE2)---All analysis will be performed on this table---------------
# * Create an object containing zeros -------------------
zeroes <- which(colSums(EcMRare2)==0)
head(zeroes)

# * Remove the columns containing zero's from OTU table ----------------------------
EcMRare3 <- EcMRare2[ , -zeroes]
head(EcMRare3)
dim(EcMRare3)#192x625

# RAREFACTION CURVES (SPECIES ACCUMULATION CURVES) ----------------------------------
SppAcumEcM <- specaccum(EcMRare3, method="exact")
plot(SppAcumEcM)




# * * EXPORT RAREFACTION TABLES  -----------------------------------------------------------------
dir.create("Analysis/QualityControl/Ectomycorrhizal/Graphs/")

#Export Rarefied table and trans2 table
write.csv(EcMRare3, "Analysis/ASV-Tables/Ectomycorrhizal/RareEcM.csv")#table w/o zeroes
write.csv(EcMTrans2, "Analysis/ASV-Tables/Ectomycorrhizal/EcMTrans2.csv")#table w/o zeroes
write.csv(EcMTotSeq, "Analysis/QualityControl/Ectomycorrhizal/EcMRepSeq.csv")


pdf("Analysis/QualityControl/Ectomycorrhizal/Graphs/MantelTest-Rarefaction.pdf", height=6, width=5)
MantelCorr
dev.off()

# Expor Species accumulation curves for comparison-------------------------------------
pdf('Analysis/QualityControl/Ectomycorrhizal/Graphs/SpeciesAccumulation.pdf', width=6, height=8)
plot(SppAcumEcM ,xlab="number of samples",ylab="number of ASV's", main="EcMRare3")
dev.off()



#***********************************************************************************************************----
#------- RAREFY THE METADATA TO MATCH RARE TABLE-------------------------------------------------------------
#***********************************************************************************************************----
## Remove samples that were removed during rarefaction from metadata =
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3

metadata<-read.csv("Metadata-PIPO.csv", na.strings = "N/A", header = TRUE) #Full Metadata

# * * LOAD OTU TABLES -----------------------------------------------------
EcMRare3<-read.csv("Analysis/ASV-Tables/Ectomycorrhizal/RareEcM.csv", check.names = FALSE) #Raw OTU table W/O Zeros (ADDED column name renamed to SampleID)
dim(EcMRare3)#192X622

#remove unwanted samples to rarefy matadata
MetaRareEcM<-metadata %>% semi_join(EcMRare3, by = "SampleID") # keep rows with matching ID
dim(MetaRareEcM)#192x26 (verity same rows as above)

#Export rarefied metadata ---------------------------------------
dir.create("Analysis/Metadata/Ectomycorrhizal/")

write.csv(MetaRareEcM, "Analysis/Metadata/Ectomycorrhizal/MetaRareEcM.csv")




##
#**************************************************************************************************************************----
# CALCULATE ALPHADIVERITY METRICS  --------------------------------------------------------------------------------------------
#**************************************************************************************************************************----

# RICHNESS METHOD 1: BIODIVERSITY CALCULATED------------------------------------------------------------------------------
# * Load metadata that correspands to R-rarefied table created above -------
attach(MetaRareEcM)
EcMRare3<-read.csv("Analysis/ASV-Tables/Ectomycorrhizal/RareEcM.csv", row.names = 1, check.names = FALSE)#reload but put SampleID as row.names


# * Richness calculations using Biodivesity package-----------------------------------
EcMRichness <- estimateR(EcMRare3)
EcMRichness<- t(estimateR(EcMRare3))#run function on transformed table to switch rows and columns

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
shanEntro <- diversity(EcMRare3) # Shannon entropy
shannon <- exp(EcMRare3 ) ## Shannon number of diversity
simpson <- diversity(EcMRare3, "inv") ## Simpson Inverse Diversity

# * Dataframe of shannon entropy, diversity &  simpson diversity--------------
EcM.richness <- data.frame(shanEntro, shannon, simpson)

# * Add above diversity metrics to S obs, chao1, ACE
EcMSppRichness <- cbind(EcM.richness, EcMRichness)
dim(EcMSppRichness)#192X628


#****************************************************************************************************************----
# EXPORT DIVERSITY DATA AND ATTACH TO RAREFIED METADATA ---------------------------------------------------
#****************************************************************************************************************----

# * * * Export Biodiversity package spp richness--------------------------------------------------------------------
dir.create("Analysis/SpeciesRichness/Ectomycorrhizal")
write.csv(EcMSppRichness, "Analysis/SpeciesRichness/Ectomycorrhizal/EcMCoreMetricsAll.csv") #all above species calculations


#sUBSET ONLY THE RICHNESS METRICS, they are at end of file and one in 1st column.................
ncol(EcMSppRichness)#628

#1 column and 6th column
CoreMetrics1<-EcMSppRichness[,623:628];head(CoreMetrics1)#628-5
CoreMetrics2<-EcMSppRichness[,1];head(CoreMetrics2)

CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)
names(CoreMetrics)[7]<- "ShannonEntropy"; head(CoreMetrics)


#Export core metrics results----------------------------------------------
dir.create("Analysis/Diversity/Ectomycorrhizal/")

write.csv(CoreMetrics, "Analysis/Diversity/Ectomycorrhizal/CoreMetrics.csv")

##
#****************************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare ---------------------------------------------------
#****************************************************************************************************************----

CoreMetrics<-read.csv("Analysis/Diversity/Ectomycorrhizal//CoreMetrics.csv", check.names = FALSE)
names(CoreMetrics)[1]<- "SampleID"


MetaRareSpp<-merge(MetaRareEcM,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#192x33

#Export rarefied metadata w core metrics attached------------------------------------------
write.csv(MetaRareSpp, "Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.csv")

