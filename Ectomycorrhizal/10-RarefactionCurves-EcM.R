#September 14, 2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
#setwd("C:/Users/fabipc/Dropbox/6-PIPO")
setwd("C:/Users/juchu/Dropbox/6-PIPO")




#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)

#LOAD QIIME DATA-NON-NORMALIZED BECUASE I WILL RAREFY IN PHYLOSEQ -------------------------------------------------------------------------------------------
Metadata<-read_tsv("Metadata/Metadata-PIPO-NC.tsv")#
Table<-read_qza("Qiime/Ectomycorrhizal/DADA2/EcM-Table-NS.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Ectomycorrhizal/Phylo/EcM-EcM-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
              c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT-----------------------------------------------------------------------------
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



#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ...........................................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor


physeqBurn<-subset_samples(physeq,Treatment=="Burned");physeqBurn #199x9116 taxa
sample_names(physeqBurn) #199 x33

physeqUnburn<-subset_samples(physeq,Treatment=="Unburned");physeqUnburn #101x9116taxa
sample_names(physeqUnburn) #101 x 33


#********************************************************************************************************************----
#------------------------------------     ESTIMATE SPECIES RICHNESS       ---------------------------------------------
#********************************************************************************************************************----
# using GlobalPatterns
library(ranacapa)
RarefyPlotEcM <- ggrare(physeq, step = 118, color = "Site",  se = TRUE)

TrtEcM <- RarefyPlotEcM + facet_wrap(~Treatment, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text( size = 16 ),
        axis.text.y = element_text( size = 16 ),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TrtEcM



SiteEcM <- RarefyPlotEcM + facet_wrap(~Site, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text( size = 16 ),
        axis.text.y = element_text( size = 16 ),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));SiteEcM

TSFplotEcM <- RarefyPlotEcM + facet_wrap(~TimeSinceFire, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text( size = 16 ),
        axis.text.y = element_text( size = 16 ),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TSFplotEcM



#*********************************************************************************************----
#--EXPORT GRAPS ----------------------------------------------------------------------------------
#*********************************************************************************************----

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/RarefactionCurves-Treatment.pdf", height=6, width=10.5)
TrtEcM
dev.off()


pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/RarefactionCurves-Site.pdf", height=10, width=10)
SiteEcM
dev.off()

pdf("Analysis/Diversity/Ectomycorrhizal/Graphs/RarefactionCurves-TSF.pdf", height=10, width=10)
TSFplotEcM
dev.off()









