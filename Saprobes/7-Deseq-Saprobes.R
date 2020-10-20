#Date last time code was ran Sept 15, 2020 

#Reset R's Ev
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO/")


#Load librarires--------------------------------------------------------------------------------------------
library("qiime2R")
library(DESeq2)#log2fold analysis
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(ochRe)# color paletteibrary(multcomp)#multiple comparisons--for glm 
library(tidyverse)



#LOAD QIIME DATA-NON-NORMALIZED ----------------------------------------------------------------------------
Metadata<-read_tsv("Metadata/Metadata-PIPO-NC.tsv")
Table<-read_qza("Qiime/Saprobes/DADA2/Saprobe-Table-NS.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Saprobes/Phylo/Saprobe-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                  c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT----------------------------------------------------------------------------------
physeqSap<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID"))
)


#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeqSap)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeqSap))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeqSap)[, "Kingdom"] <- gsub("D_0__", "", tax_table(physeqSap)[, "Kingdom"])
tax_table(physeqSap)[, "Phylum"] <- gsub("D_1__", "", tax_table(physeqSap)[, "Phylum"])
tax_table(physeqSap)[, "Class"] <- gsub("D_2__", "", tax_table(physeqSap)[, "Class"])
tax_table(physeqSap)[, "Order"] <- gsub("D_3__", "", tax_table(physeqSap)[, "Order"])
tax_table(physeqSap)[, "Family"] <- gsub("D_4__", "", tax_table(physeqSap)[, "Family"])
tax_table(physeqSap)[, "Genus"] <- gsub("D_5__", "", tax_table(physeqSap)[, "Genus"])
tax_table(physeqSap)[, "Species"] <- gsub("D_6__", "", tax_table(physeqSap)[, "Species"])


sample_data(physeqSap)$TimeSinceFire<-factor(sample_data(physeqSap)$TimeSinceFire)

#Remane taxa so it fits graph..........................................................................................
physeqSap$Genus[physeqSap$Genus == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia-CP" 

taxa_names(physeqSap)

#
#
#************************************************************************************************************----
#----------------- RUN DESEQ ANALYSIS----------------------------------------------------------------------------
#************************************************************************************************************----

#DESEQ2--------------- burned and unburned treAtment------
head(sample_data(physeqSap)$Treatment, 100)

#Ran in case there are any n/a in data
physeqSap<-subset_samples(physeqSap, Treatment!= "None")

#Import phyloseq data to create a deseq object
dds<-phyloseq_to_deseq2(physeqSap, ~ Treatment)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)

dds<-estimateSizeFactors(dds, geoMeans = geoMeans)

#Set reference level------------------------------------
dds$Treatment<-relevel(dds$Treatment,ref="Unburned")
levels(dds$Treatment)


dds<-DESeq(dds, fitType="local")


#*************************************************************************************************************----
#RESULTS----------------------------------------------------------------------------------------------------------
#*************************************************************************************************************----
res<-results(dds)
res<-res[order(res$padj, na.last=NA), ]
alpha<-0.0001 #set a cutoff value 
sigtab<-res[(res$padj < alpha), ] #create a new object containing the taxa w alpha under 0.05
sigtab<-cbind(as(sigtab, "data.frame"), as(tax_table(physeqSap)[rownames(sigtab), ], "matrix"))
sigtab


#LOOK AT TABLES MOST SIGNIFICANT FOR SITE, ONLY THE TOP ONES THAT HAD A FOLD OF 2LOGS
posigtab<-sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab<-posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
                       "Phylum", "Class", "Family", "Genus","Species")]


#
#
#************************************************************************************************************----
#Export tables---------------------------------------------------------------------------------------------------
#************************************************************************************************************----
dir.create(file.path("Analysis/Deseq/Saprobes/Tables/"), recursive = TRUE)
dir.create(file.path("Analysis/Deseq/Saprobes/Graphs/"), recursive = TRUE)

#Tables--------------------------------------------------------------------------------------
write.csv(sigtab,file="Analysis/Deseq/Saprobes/Tables/Deseq-Significant0.0001-new.csv")
write.csv(posigtab,file="Analysis/Deseq/Saprobes/Tables/Deseq-Log2Fold1.csv")


#Disperssion graph----------------
pdf("Analysis/Deseq/Saprobes/Graphs/EcM-Trt-Disperssion.pdf", height=6, width=8)
plotDispEsts(dds)
dev.off()


#
#
#************************************************************************************************************----
#RE-IMPORT REBLASTED SIGTAB FILE  AND CREATE GRAPHS--------------------------------------------------------------
#************************************************************************************************************----
#Edited file containes the same information as original file with the exception that 
#taxa that were considered unidentified were hand blasted to see if we could get a higher 
#resolution


sigtab<-read.csv("Analysis/Deseq/Saprobes/Tables/Deseq-Significant0.0001_Editted-all-info.csv", row.names=1 )#full metadata



#CREATE GRAPHING PERIMETER..........................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))



#cREATE THE GRAPH ...........................................................................
PhylumGenSap<-ggplot(sigtabgen, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
               geom_point(size=5) + 
               scale_color_manual(values=c("#45778f","#802f31","#d49b00","#c2ad97"))+ 
               theme(axis.text.x = element_text(angle = -70,hjust=0.06, 
                      vjust=.8,size=16, color = "black",face = "italic"),
                    panel.border = element_rect(color  = "black", fill=NA, size=.7), 
                    text = element_text(size=18));PhylumGenSap
            




#GRAPH 2 ---------------------------------------------------------------------------------------
#Phylum-species graph..................................................................

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen = subset(sigtab, !is.na(Species))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels = names(x))


#cREATE THE GRAPH ...........................................................................
PhylumSppSap<-ggplot(sigtabgen, aes(x = Species, y = log2FoldChange, color = Phylum)) + 
  geom_point(size=5) + 
  scale_color_manual(values=c("#45778f","#802f31","#d49b00","#c2ad97"))+ 
  theme(axis.text.x = element_text(angle = -70,hjust=0.06, 
                                   vjust=.8,size=14, color = "black",face = "italic"),
        axis.text.y = element_text(size=16, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=18));PhylumSppSap



#GRAPH 3 ---------------------------------------------------------------------------------------
#Phylum-genus graph..................................................................

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))

# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels = names(x))



#CREATE GRAPH ---------------------------------------------------------------------------------

GenSpp<-ggplot(sigtabgen, aes(x = Species, y = log2FoldChange, color = Genus)) +
  geom_point(size=5) + 
  scale_color_ochre(values=olsen_seq)+ 
  theme(axis.text.x = element_text(angle = -70,hjust=0.06, vjust=.8,size=14),
        axis.text.y = element_text(size=16),
        legend.text=element_text(size=14,face = "italic"), 
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=20));GenSpp

GenSpecies<-GenSpp+guides(col = guide_legend(nrow = 21));GenSpecies



#
#
#*************************************************************************************************----
#----- EXPORT GRAPS ----------------------------------------------------------------------------------
#*************************************************************************************************----
#
pdf("Analysis/Deseq/Saprobes/Graphs/Phylum-Genus.pdf", height=6, width=12)
PhylumGenSap
dev.off()

pdf("Analysis/Deseq/Saprobes/Graphs/Phylum-Spp.pdf", height=6, width=12)
PhylumSppSap
dev.off()

pdf("Analysis/Deseq/Saprobes/Graphs/Genus-Spp.pdf", height=8, width=12)
GenSpecies
dev.off()


#Keep file and graphs on environment as need to run the EcM codes and then create panels
# for publication







































