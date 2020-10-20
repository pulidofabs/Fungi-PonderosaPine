#Sept 15 2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO/")


#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
#library(readr)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)
library(RAM)

#LOAD QIIME DATA- RAREFIED ----------------------------------------------------------------------------
Metadata<-read_tsv("Analysis/Metadata/Saprobes/MetaRareCoreSaprobe.tsv")#Rare metadata
Table<-read_qza("Qiime/Saprobes/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Saprobes/Phylo/Saprobe-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
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
physeqBurnSap<-subset_samples(physeqSap,Treatment=="Burned")
sample_names(physeqBurnSap);physeqBurnSap #106 total

physeqUnburnSap<-subset_samples(physeqSap,Treatment=="Unburned")
sample_names(physeqUnburnSap);physeqUnburnSap #106 total


#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----

#Treatment........................................................................
Genus <- tax_glom(physeqSap, taxrank = 'Genus') # agglomerate taxa
(GenTrt = merge_samples(Genus, "Treatment")) # merge samples on sample variable of interest
RelGenTrt <- transform_sample_counts(GenTrt, function(x) x/sum(x)) #get abundance in %
RelGenTrt1sap <- psmelt(RelGenTrt) # create dataframe from phyloseq object
RelGenTrt1sap$Genus <- as.character(RelGenTrt1sap$Genus) #convert to character
RelGenTrt1sap$Genus[RelGenTrt1sap$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance

#Burned ...........................................................................
GenusB <- tax_glom(physeqBurnSap, taxrank = 'Genus') # agglomerate taxa
(GenTsfB = merge_samples(GenusB, "TimeSinceFire")) # merge samples on sample variable of interest
RelGenTsfB <- transform_sample_counts(GenTsfB, function(x) x/sum(x)) #get abundance in %
RelGenTsf1B <- psmelt(RelGenTsfB) # create dataframe from phyloseq object
RelGenTsf1B$Genus <- as.character(RelGenTsf1B$Genus) #convert to character
RelGenTsf1B$Genus[RelGenTsf1B$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance
RelGenTsf1B$Sample<-as.factor(RelGenTsf1B$Sample)
str(RelGenTsf1B$Sample)

#Unburned .......................................................................................................
GenusUn <- tax_glom(physeqUnburnSap, taxrank = 'Genus') # agglomerate taxa
(GenTsfUn = merge_samples(GenusUn, "TimeSinceFire")) # merge samples on sample variable of interest
RelGenTsfUn <- transform_sample_counts(GenTsfUn, function(x) x/sum(x)) #get abundance in %
RelGenTsf1Un <- psmelt(RelGenTsfUn) # create dataframe from phyloseq object
RelGenTsf1Un$Genus <- as.character(RelGenTsf1Un$Genus) #convert to character
RelGenTsf1Un$Genus[RelGenTsf1Un$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance



#EXPORT RELATIVE ABUNDANCE FILES .........................................................................................
dir.create(file.path("Analysis/RelativeAbundance/Saprobes/Tables/"), recursive=TRUE)

write.csv(RelGenTrt1sap, file="Analysis/RelativeAbundance/Saprobes/Tables/RelGenTrtSap-RelAbundance-2per.csv")
write.csv(RelGenTsf1B, file="Analysis/RelativeAbundance/Saprobes/Tables/RelGenTSFSap-RelAbundance-Burn-2per.csv")
write.csv(RelGenTsf1Un, file="Analysis/RelativeAbundance/Saprobes/Tables/RelGenTSFSap-RelAbundance-Unburn-2per.csv")


#
##
#********************************************************************************************************************************----
# ---------------------------------------           GRAPHS              ------------------------------------------------------
#********************************************************************************************************************************----

#Check length of colors needed............................................
length(unique(RelGenTrt1sap$Genus))#14
length(unique(RelGenTsf1B$Genus))#25
length(unique(RelGenTsf1Un$Genus))#15

#CREATE COLOR PALETTE FOR PLOTS----------------------------------------------------------------------------
#Now we generate a vector with all variable names and a vector of associated colors

All<-c("#9a9696","#3d271c","#695b4e","#986645","#ad8d5f","#a5a084","#d6cda7",
       "#ceb4a4","#7b5250","#7e716b","#9fabb7","#657b8f","#476f84","#2d5471",
       "#abad71","#8ea594","#647d53","#455849","#3e7547","#c27721","#352e4d",
       "#4b435e","#6e687e","#c4bfca","#e3d675","#a19637")

TrtSap<-c("#9a9696","#695b4e","#986645","#ad8d5f","#a5a084","#7b5250","#657b8f",
          "#abad71","#455849","#3e7547","#352e4d","#6e687e","#c4bfca","#a19637")

Burn<-c("#9a9696","#3d271c","#695b4e","#986645","#ad8d5f","#d6cda7","#ceb4a4",
        "#7b5250","#7e716b","#657b8f","#476f84","#24455e", "#abad71","#8ea594",
        "#647d53","#455849","#3e7547","#c4a700", "#c27721","#352e4d","#4b435e",
        "#6e687e","#c4bfca","#e3d675","#a19637")

Un<-c("#9a9696","#695b4e","#a5a084","#d6cda7","#ceb4a4","#7b5250","#9fabb7","#657b8f",
      "#abad71","#455849","#3e7547","#6e687e","#c4bfca","#a19637","#373c40")


#TREATMENT-------------------------------------------------------- # NOT SURE Y THEY ARE AT 200------------
AbunTrtSap<-ggplot(data=RelGenTrt1sap, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Burn)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=18), 
                    axis.text.y = element_text(size=18, hjust=0.5),
                    axis.text.x = element_text(size=18), 
                    legend.position="right",
                    legend.text = element_text(face = "italic",size = 14),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  xlab("Treatment")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol=1));AbunTrtSap


#BURNED PLOTS-----------------------------------------------------------------------------------
AbunSapB<-ggplot(data=RelGenTsf1B, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Burn)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=20), 
                    axis.text.y = element_text(size=16, hjust=0.5),
                    axis.text.x = element_text(size=16),
                    legend.position="right",
                    legend.text = element_text(face = "italic", size=16),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunSapB


#BURNED PLOTS------------------------------------------------------------------------------------
AbunSapUn<-ggplot(data=RelGenTsf1Un, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Un)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=20), 
                    axis.text.y = element_text(size=16, hjust=0.5),
                    axis.text.x = element_text(size=16),
                    legend.position="right", 
                    legend.text = element_text(face = "italic", size = 16),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol=1));AbunSapUn


#Export Graphs-------------------------------------------------------------------------------------
GenTSF<-ggarrange(AbunSapB, AbunSapUn,  ncol=1, nrow=2, labels= c("A Burned","B Unburned"), 
                  font.label = list(size = 16,face = "bold", color ="black"),
                  common.legend = TRUE,legend = "right", align = "hv");GenTSF


GenTSF2<-ggarrange(AbunSapB, AbunSapUn,  ncol=2, nrow=1, labels= c("A Burned","B Unburned"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   legend = "right", align = "hv");GenTSF2





#********************************************************************************************************----
#Export Graphs-----------------------------------------------------------------------------------------------
#********************************************************************************************************----
dir.create(file.path("Analysis/RelativeAbundance/Saprobes/Graphs/"), recursive=TRUE)


pdf("Analysis/RelativeAbundance/Saprobes/Graphs/Rel-Abund-Trt.pdf", height=5, width=6.5)
AbunTrtSap
dev.off()

pdf("Analysis/RelativeAbundance/Saprobes/Graphs/Rel-Abund-TSF-Burned.pdf", height=5, width=6.5)
AbunSapB 
dev.off()

pdf("Analysis/RelativeAbundance/Saprobes/Graphs/Rel-Abund-TSF-Unburned.pdf", height=5, width=6.5)
AbunSapUn
dev.off()

#Export Panels--------------------------------------------------------------------------------------------
pdf("Analysis/RelativeAbundance/Saprobes/Graphs/TSF-Burned-Unburned1a.pdf", height=8, width=7)
GenTSF
dev.off()

pdf("Analysis/RelativeAbundance/Saprobes/Graphs/TSF-Burned-Unburned2.pdf", height=4, width=8)
GenTSF2
dev.off()


