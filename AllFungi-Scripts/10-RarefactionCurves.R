#Last time code was ran............ June 20, 2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO")


#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)

#LOAD QIIME DATA-NON-NORMALIZED -------------------------------------------------------------------------------------------
Metadata<-read_tsv("Metadata-PIPO-NC.tsv")#Rare metadata
Table<-read_qza("QiimeOutputs/DADA2/Fungal-Table-dada2-NS.qza", tmp = "C:/tmp")#Need unrarefied table
Tree<-read_qza("QiimeOutputs/Phylo/Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("QiimeOutputs/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
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
RarefyPlot <- ggrare(physeq, step = 6229, color = "Site",  se = TRUE)
             
Trt <- RarefyPlot + facet_wrap(~Treatment, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 16 ),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
       legend.position="bottom",strip.text = element_text(size = 18));Trt


Site <- RarefyPlot + facet_wrap(~Site, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
       axis.text.x = element_text(size = 15),
       axis.text.y = element_text(size = 16),
       panel.background = element_rect(fill = 'white', colour = 'black'),
       panel.spacing = unit(1, "lines"),
       legend.position="bottom",strip.text = element_text(size = 18));Site


TSFplot <- RarefyPlot + facet_wrap(~TimeSinceFire, scales = "free_x")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TSFplot



#*********************************************************************************************----
#--EXPORT GRAPS ----------------------------------------------------------------------------------
#*********************************************************************************************----

pdf("Analysis/Diversity/Fungi/Graphs/RarefactionCurves-Treatment.pdf", height=6, width=10.5)
Trt
dev.off()


pdf("Analysis/Diversity/Fungi/Graphs/RarefactionCurves-Site.pdf", height=10, width=10)
Site
dev.off()

pdf("Analysis/Diversity/Fungi/Graphs/RarefactionCurves-TSF.pdf", height=10, width=10)
TSFplot
dev.off()



#*********************************************************************************************----
#--EXPORT PANELS ----------------------------------------------------------------------------------
#*********************************************************************************************----
#ran each code individually then combined them with this code
#


pdf("Analysis/Diversity/RarefactionCurve-Trt.pdf", height=12, width=10, onefile=FALSE)
ggarrange(Trt,TrtEcM,TrtSap, ncol=1, nrow=3, common.legend = TRUE, 
          labels = c("a Fungi","b EcM", "c Saprobes"),
          legend="bottom", hjust = c(-.3, .1),
          align="hv",font.label = list(size=16,face="bold"))
dev.off()

pdf("Analysis/Diversity/RarefactionCurve-TSF-EcM-Sap.pdf", height=12, width=10, onefile=FALSE)
ggarrange(TSFplotEcM,TSFplotSap, ncol=1, nrow=2, common.legend = TRUE, 
          labels = c("a EcM", "b Saprobes"),legend="bottom", hjust = c(-.3, .1),
          align="hv",font.label = list(size=16,face="bold"))
dev.off()