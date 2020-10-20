#September 14, 2020

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

#LOAD QIIME DATA-RAREFIED----------------------------------------------------------------------------
Metadata<-read_tsv("Analysis/Metadata/Ectomycorrhizal/MetaRareCoreEcM.tsv")#
Table<-read_qza("Qiime/Ectomycorrhizal/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
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
#Quality control: Remove the g__ from each rank number-Burned
rank_names(physeq)
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])

sample_data(physeq)$TimeSinceFire<-factor(sample_data(physeq)$TimeSinceFire)



#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ------------------------------------------
physeqBurn<-subset_samples(physeq,Treatment=="Burned")
sample_names(physeqBurn);physeqBurn #106 total

physeqUnburn<-subset_samples(physeq,Treatment=="Unburned")
sample_names(physeqUnburn);physeqUnburn #106 total


#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----

#Treatment........................................................................
Genus <- tax_glom(physeq, taxrank = 'Genus') # agglomerate taxa
(GenTrt = merge_samples(Genus, "Treatment")) # merge samples on sample variable of interest
RelGenTrt <- transform_sample_counts(GenTrt, function(x) x/sum(x)) #get abundance in %
RelGenTrt1 <- psmelt(RelGenTrt) # create dataframe from phyloseq object
RelGenTrt1$Genus <- as.character(RelGenTrt1$Genus) #convert to character
RelGenTrt1$Genus[RelGenTrt1$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance


#Burned ...........................................................................
GenusB <- tax_glom(physeqBurn, taxrank = 'Genus') # agglomerate taxa
(GenTsfB = merge_samples(GenusB, "TimeSinceFire")) # merge samples on sample variable of interest
RelGenTsfB <- transform_sample_counts(GenTsfB, function(x) x/sum(x)) #get abundance in %
RelGenTsf1B <- psmelt(RelGenTsfB) # create dataframe from phyloseq object
RelGenTsf1B$Genus <- as.character(RelGenTsf1B$Genus) #convert to character
RelGenTsf1B$Genus[RelGenTsf1B$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance
RelGenTsf1B$Sample<-as.factor(RelGenTsf1B$Sample)


#Unburned ......................................................................
GenusUn <- tax_glom(physeqUnburn, taxrank = 'Genus') # agglomerate taxa
(GenTsfUn = merge_samples(GenusUn, "TimeSinceFire")) # merge samples on sample variable of interest
RelGenTsfUn <- transform_sample_counts(GenTsfUn, function(x) x/sum(x)) #get abundance in %
RelGenTsf1Un <- psmelt(RelGenTsfUn) # create dataframe from phyloseq object
RelGenTsf1Un$Genus <- as.character(RelGenTsf1Un$Genus) #convert to character
RelGenTsf1Un$Genus[RelGenTsf1Un$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance


#EXPORT RELATIVE ABUNDANCE TABLES----------------------------------------------------------------------------------------------------
write.csv(RelGenTrt1, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundaceTrt-2per.csv")
write.csv(RelGenTsf1B, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundaceTimeSinceFire-Burned-2per.csv")
write.csv(RelGenTsf1Un, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundaceTimeSinceFire-Unburned-2per.csv")




#
#********************************************************************************************************************************----
# ----------------------------------------------  GRAPHS  ---------------------------------------------------------------------------
#********************************************************************************************************************************----

#Check length of colors needed...............................................................
length(unique(RelGenTrt1$Genus))
length(unique(RelGenTsf1B$Genus))
length(unique(RelGenTsf1Un$Genus))

#CREATE COLOR PALETTE FOR PLOTS...............................................................
name.check((RelGenTrt1$Genus))
levels(RelGenTrt1$Genus)

#Now we generate a vector with all variable names and a vector of associated colors

Trt<-c("< 2% abund."="#9a9696","Cenococcum"="#470000","Chromelosporium"="#8e3600",
       "Cortinarius"="#c25715","Helvellosebacina"="#d98b5b","Hygrophorus"="#7a5334",
      "Inocybe"="#382752","Lyophyllum"="#73637d","Peziza"="#a59bab","Piloderma"="#d5cede",
      "Pustularia"="#c2bc8a","Rhizopogon"="#f2efd3","Russula"="#e1f0fd",
      "Thelephora"="#aec4ce","Tomentella"="#235494","Tuber"="#445b39",
      "Wilcoxina"="#38381c")

Burn<-c("< 2% abund."="#9a9696","Cenococcum"="#470000","Chromelosporium"="#8e3600",
        "Cortinarius"="#c25715", "Geopyxis"="#b08c6f","Hebeloma"="#d4bf3b",
        "Inocybe"="#382752","Lyophyllum"="#73637d","Morchella"="#9e6c6c",
        "Peziza"="#a59bab","Pustularia"="#c2bc8a","Rhizopogon"="#f2efd3",
        "Sistotrema"="#c7c7c7","Thelephora"="#aec4ce","Trichophaea"="#466343",
        "Wilcoxina"="#38381c")
 

Un<-c("< 2% abund."="#9a9696","Amphinema"="#9e0022","Cortinarius"="#c25715",
      "Entoloma"="#c4835a","Helvellosebacina"="#d98b5b","Hygrophorus"="#7a5334",
      "Inocybe"="#382752","Piloderma"="#d5cede","Rhizopogon"="#f2efd3",
      "Russula"="#e1f0fd","Sistotrema"="#c7c7c7","Suillus"="#8c91cf",
      "Thelephora"="#aec4ce","Tomentella"="#235494","Tricholoma"="#84a38c",
      "Tuber"="#647d5c","Tylospora"="#056e20","Wilcoxina"="#38381c")
      
      

#TREATMENT------------------------------------------------------------------------------
AbunTrt<-ggplot(data=RelGenTrt1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Trt)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=18), 
        axis.text.y = element_text(size=18, hjust=0.5),
        axis.text.x = element_text(size=18), 
        legend.position="right",
        legend.text = element_text(face = "italic",size = 14),
        plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  xlab("Treatment")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol=1));AbunTrt
  

#BURNED PLOTS-----------------------------------------------------------------------------------
AbunB<-ggplot(data=RelGenTsf1B, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Burn)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=20), 
        axis.text.y = element_text(size=16, hjust=0.5),
        axis.text.x = element_text(size=16),
        legend.position="right",
        legend.text = element_text(face = "italic", size=16),
        plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunB
  


#BURNED PLOTS------------------------------------------------------------------------------------
AbunUn<-ggplot(data=RelGenTsf1Un, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Un)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(size=20), 
        axis.text.y = element_text(size=16, hjust=0.5),
        axis.text.x = element_text(size=16),
        legend.position="right",
        legend.text = element_text(face = "italic", size=16),
        plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunUn


#Export Graphs-------------------------------------------------------------------------------------
GenTSF<-ggarrange(AbunB, AbunUn,  ncol=1, nrow=2, labels= c("a Burned","b Unburned"), 
                  font.label = list(size = 16,face = "bold", color ="black"),
                  common.legend = FALSE,legend = "right", align = "hv");GenTSF


GenTSF2<-ggarrange(AbunB, AbunUn,  ncol=2, nrow=1, labels= c("a Burned","b Unburned"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF2


#Export Graphs--------------------------------------------------------------------------------------------
pdf("Analysis/RelativeAbundance/Ectomycorrhizal/Graphs/Rel-Abund-Trt-New.pdf", height=5, width=6.5)
AbunTrt
dev.off()

pdf("Analysis/RelativeAbundance/Ectomycorrhizal/Graphs/Rel-Abund-TSF-Burned-New.pdf", height=5, width=6.5)
AbunB
dev.off()

pdf("Analysis/RelativeAbundance/Ectomycorrhizal/Graphs/Rel-Abund-TSF-Unburned-New.pdf", height=5, width=6.5)
AbunUn
dev.off()


#Export Panels--------------------------------------------------------------------------------------------
pdf("Analysis/RelativeAbundance/Ectomycorrhizal/Graphs/TSF-Burned-Unburned-New.pdf", height=10, width=7)
GenTSF
dev.off()


pdf("Analysis/RelativeAbundance/Ectomycorrhizal/Graphs/TSF-Burned-Unburned-New.pdf", height=6, width=14)
GenTSF2
dev.off()




#PANELS OF ECM AND SAPROBE BURNED AND UNBURNED--------------------------------------------------------------
#Ran the relative abundance individually for EcM and Sap and used the script below to create the panels
#

GenTSF1<-ggarrange(AbunTrt,AbunTrtSap,ncol=2, nrow=1, labels= c("a EcM","b Saprobes"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF1


GenTSF2<-ggarrange(AbunB,AbunSapB,ncol=1, nrow=2, labels= c("a EcM","b Saprobes"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF2


GenTSF3<-ggarrange(AbunB, AbunUn,AbunSapB,AbunSapUn, ncol=2, nrow=2, 
                   labels= c("a Buned","b Unburned", "c Burned", "d Unburned"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF3



pdf("Analysis/RelativeAbundance/Panels/Rel-Abund-Treatment-New.pdf", height=10, width=6)
GenTSF1
dev.off()




pdf("Analysis/RelativeAbundance/Panels/Rel-Abund-TSF-ALL-New.pdf", height=12, width=14)
GenTSF3
dev.off()

pdf("Analysis/RelativeAbundance/Panels/Rel-Abund-TSF-Burned-New.pdf", height=12, width=7)
GenTSF2
dev.off()









