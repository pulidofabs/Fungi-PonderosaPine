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

#LOAD QIIME DATA-RAREFIED---------------------------------------------------------------------------------------------------
Metadata<-read_tsv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.tsv")#
Table<-read_qza("QiimeOutputs/Ectomycorrhizal/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
Tree<-read_qza("QiimeOutputs/Ectomycorrhizal/Phylo/EcM-EcM-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("QiimeOutputs/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                 c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT---------------------------------------------------------------------------------------------------
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID"))
)



#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])



#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ------------------------------------------
sample_data(physeq)$TimeSinceFire<-factor(sample_data(physeq)$TimeSinceFire) # make as factor

physeqBurn<-subset_samples(physeq,Treatment=="Burned");physeqBurn #87x604 taxa 
sample_names(physeqBurn) #199 x33

physeqUnburn<-subset_samples(physeq,Treatment=="Unburned");physeqUnburn #105x604
sample_names(physeqUnburn) #101 x 33



#.----
#.----
#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----

#Treatment........................................................................
Genus <- tax_glom(physeq, "Genus")
RelGenTrt <- transform_sample_counts(Genus, function(x) x / sum(x))
RelGenTrt1<-psmelt(RelGenTrt)

#Burned ...........................................................................
GenusB<- tax_glom(physeqBurn, "Genus")
RelGenTsfB <- transform_sample_counts(GenusB, function(x) x / sum(x))
RelGenTsf1B<-psmelt(RelGenTsfB)

#Unburned ......................................................................
GenusUn<- tax_glom(physeqUnburn,"Genus")
RelGenTsfUn <- transform_sample_counts(GenusUn, function(x) x / sum(x))
RelGenTsf1Un<-psmelt(RelGenTsfUn)


#.----
#.----
#**************************************************************************************************************************-----
# QUALITY CONTROL FOR SUBSET DATASET  ------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----

#Set reference level to unburn..........................................................................
RelGenTrt1$Treatment<-as.factor(RelGenTrt1$Treatment)
RelGenTrt1$Treatment <- try(relevel(RelGenTrt1$Treatment , "Unburned"));levels(RelGenTrt1$Treatment)


RelGenTsf1B$TimeSinceFire<-as.factor(RelGenTsf1B$TimeSinceFire)
RelGenTsf1B$TimeSinceFire <- try(relevel(RelGenTsf1B$TimeSinceFire , "2"));levels(RelGenTsf1B$TimeSinceFire)


RelGenTsf1Un$TimeSinceFire<-as.factor(RelGenTsf1Un$TimeSinceFire)
RelGenTsf1Un$TimeSinceFire <- try(relevel(RelGenTsf1Un$TimeSinceFire , "2"));levels(RelGenTsf1Un$TimeSinceFire)


#EXPORT FILES ....................................................................................................
dir.create("Analysis/RelativeAbundance/Ectomycorrhizal/Tables")

write.csv(RelGenTrt1, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundaceTrt-DesStats-Data.csv")
write.csv(RelGenTsf1B, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundace-TSF-Burn-DesStats-Data.csv")
write.csv(RelGenTsf1Un, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundace-TSF-Unburn-DesStats-Data.csv")




#
#
#********************************************************************************************************************----
#---DESCRIPTIVE STATISTICS -------------------------------------------------------------------------------------------------
#********************************************************************************************************************----

#Treatment..............................................................
TrtDesStats<-RelGenTrt1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment) %>%
  summarise(n_obs = n(),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TrtDesStats



#Burned.................................................................
TSFDesStatsB<-RelGenTsf1B %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TSFDesStatsB



#Unburned...............................................................
TSFDesStatsUn<-RelGenTsf1Un %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(n_obs = n(),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TSFDesStatsUn



#Export descriptive statistics........................................................................................................
dir.create("Analysis/RelativeAbundance/Ectomycorrhizal/Tables/DesStats")

write.csv(TrtDesStats, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/DesStats/TRT-DescripStats.csv")
write.csv(TSFDesStatsB, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/DesStats/TSF-Burned-DescripStats.csv")
write.csv(TSFDesStatsUn, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/DesStats/TSF-Unburned-DescripStats.csv")


#
#
#***********************************************************************************************************************----
#--- PERCENT CHANGE CALCULATIONS -------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
library(dplyr)

TrtPerChg<-RelGenTrt1%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TrtPerChg

TsfPerChgB<-RelGenTsf1B%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgB

TsfPerChgUn<-RelGenTsf1Un%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TimeSinceFire) %>%
  summarise(mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgUn



#Export tables................................................................................................
dir.create("Analysis/RelativeAbundance/Ectomycorrhizal/Tables/PercentChange")

write.csv(TrtPerChg,"Analysis/RelativeAbundance/Ectomycorrhizal/Tables/PercentChange/TRT-PerChange.csv")
write.csv(TsfPerChgB,"Analysis/RelativeAbundance/Ectomycorrhizal/Tables/PercentChange/TSF-Burned-PerChange.csv")
write.csv(TsfPerChgUn,"Analysis/RelativeAbundance/Ectomycorrhizal/Tables/PercentChange/TSF-Unburned-PerChange.csv")





