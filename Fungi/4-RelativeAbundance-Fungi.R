#9-15-2020

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

#LOAD QIIME DATA-NON-NORMALIZED -------------------------------------------------------------------------------------------
Metadata<-read_tsv("Analysis/Metadata/Fungi/MetaRareCore.tsv")#Rare metadata
Table<-read_qza("Qiime/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")#Need unrarefied table
Tree<-read_qza("Qiime/Phylo/Rooted-Tree.qza",tmp = "C:/tmp" )
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
str(RelGenTsf1B$Sample)

#Unburned ......................................................................
GenusUn <- tax_glom(physeqUnburn, taxrank = 'Genus') # agglomerate taxa
(GenTsfUn = merge_samples(GenusUn, "TimeSinceFire")) # merge samples on sample variable of interest
RelGenTsfUn <- transform_sample_counts(GenTsfUn, function(x) x/sum(x)) #get abundance in %
RelGenTsf1Un <- psmelt(RelGenTsfUn) # create dataframe from phyloseq object
RelGenTsf1Un$Genus <- as.character(RelGenTsf1Un$Genus) #convert to character
RelGenTsf1Un$Genus[RelGenTsf1Un$Abundance < 0.02] <- "< 2% abund." #rename genera with < 1% abundance


#EXPORT RELATIVE ABUNDANCE TABLES----------------------------------------------------------------------------------------------------
write.csv(RelGenTrt1, "Analysis/RelativeAbundance/Fungi/Tables/RelativeAbundaceTrt-2per.csv")
write.csv(RelGenTsf1B, "Analysis/RelativeAbundance/Fungi/Tables/RelativeAbundaceTimeSinceFire-Burned-2per.csv")
write.csv(RelGenTsf1Un, "Analysis/RelativeAbundance/Fungi/Tables/RelativeAbundaceTimeSinceFire-Unburned-2per.csv")




#********************************************************************************************************************************----
# ---------------------------------------           GRAPHS              ------------------------------------------------------
#********************************************************************************************************************************----

#Check length of colors needed............................................
length(unique(RelGenTrt1$Genus))#19
length(unique(RelGenTsf1B$Genus))#24
length(unique(RelGenTsf1Un$Genus))#24

#CREATE COLOR PALETTE FOR PLOTS----------------------------------------------------------------------------
#Now we generate a vector with all variable names and a vector of associated colors

Trt<-c("#9a9696","#4f4339","#695b4e","#986645","#ad8d5f","#c25715","#7b5250","#657b8f",
       "#e6ba0b","#bf9bae","#455849","#71637a","#352e4d","#dceded","#c4bfca","#a19637",
       "#006691","#38381c","#38381c")
       
Burn<-c("#9a9696","#4f4339","#695b4e","#986645","#ad8d5f","#d6cda7","#7b5250",
        "#657b8f","#476f84","#2d5471","#bf6f6b","#8ea594","#647d53","#58992c",
        "#657b8f","#3e7547","#d4d4d4","#574b80","#6e687e","#c4bfca","#e3d675",
        "#a19637","#006691","#38381c")

Unburn<-c("#9a9696","#695b4e","#c25715","#8f6a4f","#7b5250","#d477a8","#657b8f",
          "#c78b1b","#e6ba0b","#bf9bae","#455849","#71637a","#3e7547","#f2efd3",
          "#dceded","#6e687e","#abc6c9","#c4bfca","#7491b4","#486486","#235494",
          "#647d53","#a19637","#38381c")


#TREATMENT-------------------------------------------------------- # NOT SURE Y THEY ARE AT 200------------
AbunTrt<-ggplot(data=RelGenTrt1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Trt)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=18), 
                    axis.text.y = element_text(size=18, hjust=0.5),
                    axis.text.x = element_text(size=18), 
                    legend.position="right",
                    legend.text = element_text(face = "italic",size = 14),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  xlab("Treatment")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol=1));AbunTrt
  

#BURNED PLOTS-----------------------------------------------------------------------------------
AbunB<-ggplot(data=RelGenTsf1B, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Burn)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=20), 
                    axis.text.y = element_text(size=16, hjust=0.5),
                    axis.text.x = element_text(size=16),
                    legend.position="right",
                    legend.text = element_text(face = "italic", size=12.5),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunB
  


#BURNED PLOTS------------------------------------------------------------------------------------
AbunUn<-ggplot(data=RelGenTsf1Un, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Unburn)+
  theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size=20), 
                    axis.text.y = element_text(size=16, hjust=0.5),
                    axis.text.x = element_text(size=16),
                    legend.position="right", 
                    legend.text = element_text(face = "italic", size = 12.5),
                    plot.margin = unit(c(2,.5,1.5,0.2), "lines"))+ 
  scale_x_discrete(limits=c("2","3","5","11"))+
  xlab("Time Since Fire (Yrs)")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunUn


#Export Graphs-------------------------------------------------------------------------------------
GenTSF<-ggarrange(AbunB, AbunUn,  ncol=1, nrow=2, labels= c("A Burned","B Unburned"), 
                  font.label = list(size = 16,face = "bold", color ="black"),
                  common.legend = FALSE,legend = "right", align = "hv");GenTSF


GenTSF2<-ggarrange(AbunB, AbunUn,  ncol=2, nrow=1, labels= c("A Burned","B Unburned"), 
                   font.label = list(size = 16,face = "bold", color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF




#Export Graphs--------------------------------------------------------------------------------------------
dir.create(file.path("Analysis/RelativeAbundance/Fungi/Graphs/"), recursive = TRUE)

pdf("Analysis/RelativeAbundance/Fungi/Graphs/Rel-Abund-Trt.pdf", height=6, width=6)
AbunTrt
dev.off()

pdf("Analysis/RelativeAbundance/Fungi/Graphs/Rel-Abund-TSF-Burned.pdf", height=7.5, width=6.5)
AbunB
dev.off()

pdf("Analysis/RelativeAbundance/Fungi/Graphs/Rel-Abund-TSF-Unburned.pdf", height=7.5, width=6.5)
AbunUn
dev.off()

#Export Panels--------------------------------------------------------------------------------------------
pdf("Analysis/RelativeAbundance/Fungi/Graphs/TSF-Burned-Unburned.pdf", height=13, width=8)
GenTSF
dev.off()


pdf("Analysis/RelativeAbundance/Fungi/Graphs/TSF-Burned-Unburned-2.pdf", height=8, width=12)
GenTSF2
dev.off()










