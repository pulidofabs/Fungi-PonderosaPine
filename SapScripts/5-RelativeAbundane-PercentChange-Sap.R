#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO/")

#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)

#LOAD QIIME DATA- --------------------------------------------------------------------------------------------------------
Metadata<-read_tsv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.tsv")#Rare metadata
Table<-read_qza("QiimeOutputs/Saprobes/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
Tree<-read_qza("QiimeOutputs/Saprobes/Phylo/Saprobe-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("QiimeOutputs/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT-----------------------------------------------------------------------------
physeqSap<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID"))
)



#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

#Quality control: Remove the g__ from each rank number-Burned

rank_names(physeqSap)
colnames(tax_table(physeqSap))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeqSap)[, "Phylum"] <- gsub("p__", "", tax_table(physeqSap)[, "Phylum"])
tax_table(physeqSap)[, "Class"] <- gsub("c__", "", tax_table(physeqSap)[, "Class"])
tax_table(physeqSap)[, "Order"] <- gsub("o__", "", tax_table(physeqSap)[, "Order"])
tax_table(physeqSap)[, "Family"] <- gsub("f__", "", tax_table(physeqSap)[, "Family"])
tax_table(physeqSap)[, "Genus"] <- gsub("g__", "", tax_table(physeqSap)[, "Genus"])
tax_table(physeqSap)[, "Species"] <- gsub("s__", "", tax_table(physeqSap)[, "Species"])

sample_data(physeqSap)$TimeSinceFire<-factor(sample_data(physeqSap)$TimeSinceFire)



#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ------------------------------------------
sample_data(physeqSap)$TimeSinceFire<-factor(sample_data(physeqSap)$TimeSinceFire) # make as factor

physeqBurnSap<-subset_samples(physeqSap,Treatment=="Burned");physeqBurnSap #106x2583 taxa 
sample_names(physeqBurnSap) #199 x33

physeqUnburnSap<-subset_samples(physeqSap,Treatment=="Unburned");physeqUnburnSap #102x2583
sample_names(physeqUnburnSap) #101 x 33



#.----
#.----
#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----
#Treatment........................................................................
GenusSap <- tax_glom(physeqSap, "Genus")
RelGenTrtSap <- transform_sample_counts(GenusSap, function(x) x / sum(x))
RelGenTrt1sap<-psmelt(RelGenTrtSap)

#Burned ...........................................................................
GenusBsap<- tax_glom(physeqBurnSap, "Genus")
RelGenTsfBsap <- transform_sample_counts(GenusBsap, function(x) x / sum(x))
RelGenTsf1Bsap<-psmelt(RelGenTsfBsap)

#Unburned ........................................................................
GenusUnSap<- tax_glom(physeqUnburnSap,"Genus")
RelGenTsfUnSap <- transform_sample_counts(GenusUnSap, function(x) x / sum(x))
RelGenTsf1UnSap<-psmelt(RelGenTsfUnSap)




#.----
#.----
#**************************************************************************************************************************-----
# QUALITY CONTROL FOR SUBSET DATASET  ------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----

#Set reference level to unburn..........................................................................
RelGenTrt1sap$Treatment<-as.factor(RelGenTrt1sap$Treatment)
RelGenTrt1sap$Treatment <- try(relevel(RelGenTrt1sap$Treatment , "Unburned"));levels(RelGenTrt1sap$Treatment)


RelGenTsf1Bsap$TimeSinceFire<-as.factor(RelGenTsf1Bsap$TimeSinceFire)
RelGenTsf1Bsap$TimeSinceFire <- try(relevel(RelGenTsf1Bsap$TimeSinceFire , "2"));levels(RelGenTsf1Bsap$TimeSinceFire)


RelGenTsf1UnSap$TimeSinceFire<-as.factor(RelGenTsf1UnSap$TimeSinceFire)
RelGenTsf1UnSap$TimeSinceFire <- try(relevel(RelGenTsf1UnSap$TimeSinceFire , "2"));levels(RelGenTsf1UnSap$TimeSinceFire)


#EXPORT FILES ....................................................................................................
dir.create("Analysis/RelativeAbundance/Saprobes/Tables")

write.csv(RelGenTrt1sap, "Analysis/RelativeAbundance/Saprobes/Tables/RelativeAbundaceTrt-DesStats-Data.csv")
write.csv(RelGenTsf1Bsap, "Analysis/RelativeAbundance/Saprobes/Tables/RelativeAbundace-TSF-Burn-DesStats-Data.csv")
write.csv(RelGenTsf1UnSap, "Analysis/RelativeAbundance/Saprobes/Tables/RelativeAbundace-TSF-Unburn-DesStats-Data.csv")




#
#
#********************************************************************************************************************----
#---DESCRIPTIVE STATISTICS -------------------------------------------------------------------------------------------------
#********************************************************************************************************************----

#Treatment..............................................................
TrtDesStatsSap<-RelGenTrt1sap %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment) %>%
  summarise(n_obs = n(),
            Total = sum(Abundance, na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TrtDesStatsSap



#Burned.................................................................
TSFDesStatsBsap<-RelGenTsf1Bsap %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(n_obs = n(),
            Total = sum(Abundance, na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TSFDesStatsBsap



#Unburned...............................................................
TSFDesStatsUnSap<-RelGenTsf1UnSap %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(n_obs = n(),
            Total = sum(Abundance, na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TSFDesStatsUnSap



#Export descriptive statistics........................................................................................................
dir.create("Analysis/RelativeAbundance/Saprobes/Tables/DesStats")

write.csv(TrtDesStatsSap, "Analysis/RelativeAbundance/Saprobes/Tables/DesStats/TRT-DescripStats.csv")
write.csv(TSFDesStatsBsap, "Analysis/RelativeAbundance/Saprobes/Tables/DesStats/TSF-Burned-DescripStats.csv")
write.csv(TSFDesStatsUnSap, "Analysis/RelativeAbundance/Saprobes/Tables/DesStats/TSF-Unburned-DescripStats.csv")


#
#
#***********************************************************************************************************************----
#--- PERCENT CHANGE CALCULATIONS -------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
library(dplyr)

TrtPerChgSap<-RelGenTrt1sap%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TrtPerChgSap

TsfPerChgBsap<-RelGenTsf1Bsap%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgBsap

TsfPerChgUnSap<-RelGenTsf1UnSap%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgUnSap



#Export tables................................................................................................
dir.create("Analysis/RelativeAbundance/Saprobes/Tables/PercentChange")

write.csv(TrtPerChgSap,"Analysis/RelativeAbundance/Saprobes/Tables/PercentChange/TRT-PerChange.csv")
write.csv(TsfPerChgBsap,"Analysis/RelativeAbundance/Saprobes/Tables/PercentChange/TSF-Burned-PerChange.csv")
write.csv(TsfPerChgUnSap,"Analysis/RelativeAbundance/Saprobes/Tables/PercentChange/TSF-Unburned-PerChange.csv")





