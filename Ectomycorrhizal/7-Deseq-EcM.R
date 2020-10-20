#Date code was last ran------June 5, 2020

#Reset R'----------------
rm(list=ls())

# SET WORKING DIRECTORY-----------------------------------------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/6-PIPO")

#Load Packages--------------------------------------------------------------------------------------------------
library("qiime2R")
library(DESeq2)#log2fold analysis
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(ochRe)# color paletteibrary(multcomp)#multiple comparisons--for glm 
library(tidyverse)


#LOAD QIIME DATA-NON-NORMALIZED --------------------------------------------------------------------------------
metadata<-read_tsv("Metadata/Metadata-PIPO-NC.tsv")#full metadata
Table<-read_qza("Qiime/Ectomycorrhizal/DADA2/EcM-Table-NS.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Ectomycorrhizal/Phylo/EcM-EcM-Rooted-Tree.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
           c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT---------------------------------------------------------------------------------------
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(metadata%>% as.data.frame() %>% column_to_rownames("SampleID"))
)



#
#
#***********************************************************************************************************----
#QUALITY CONTROL------------------------------------------------------------------------------------------------
#***********************************************************************************************************----

#Quality control: Remove the Preceding letters from ranks
rank_names(physeq)

tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])


#
#
#************************************************************************************************************----
#----------------- RUN DESEQ ANALYSIS----------------------------------------------------------------------------
#************************************************************************************************************----

#DESEQ2--------------- burned and unburned treAtment------
head(sample_data(physeq)$Treatment, 100)

#Ran in case there are any n/a in data
Physeq<-subset_samples(physeq, Treatment!= "None")


#Import phyloseq data to create a deseq object
dds<-phyloseq_to_deseq2(physeq, ~ Treatment)

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
alpha<-0.05 #set a cutoff value 
sigtab<-res[(res$padj < alpha), ] #create a new object containing the taxa w alpha under 0.05
sigtab<-cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
sigtab


#LOOK AT TABLES MOST SIGNIFICANT FOR SITE, ONLY THE TOP ONES THAT HAD A FOLD OF 2LOGS
posigtab<-sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab<-posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", 
                       "Family", "Genus","Species")]


#************************************************************************************************************----
#Export tables---------------------------------------------------------------------------------------------------
#************************************************************************************************************----

#Tables-......................................................................................
dir.create(file.path("Analysis/Deseq/Ectomycorrhizal/Tables"), recursive = TRUE)
dir.create(file.path("Analysis/Deseq/Ectomycorrhizal/Graphs"), recursive = TRUE)

write.csv(sigtab,file="Analysis/Deseq/Ectomycorrhizal/Tables/Deseq-Significant0.05.csv")
write.csv(posigtab,file="Analysis/Deseq/Ectomycorrhizal/Tables/Deseq-Log2Fold.csv")


#Disperssion graph..................
pdf("Analysis/Deseq/Ectomycorrhizal/Graphs/EcM-Trt-Disperssion.pdf", height=6, width=8)
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

sigtab<-read.csv("Analysis/Deseq/Ectomycorrhizal/Tables/Deseq-Significant0.05-Editted.csv", row.names=1)



#CREATE GRAPHING PERIMETER...............................................................
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



#cREATE THE GRAPH ........................................................................
PhylumGen<-ggplot(sigtabgen, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
    geom_point(size=5) + 
    scale_color_manual(values=c("#45778f","#802f31"))+ 
    theme(axis.text.x = element_text(angle = -70,hjust=0.06, 
          vjust=.8,size=16, color = "black",face = "italic"),
          panel.border = element_rect(color  = "black", fill=NA, size=.7), 
          text = element_text(size=18));PhylumGen
   

#GRAPH 2 -----------------------------------------------------------------------------------------

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



#Plot graph..............................................................................
PhylumSpp<-ggplot(sigtabgen, aes(x = Species, y = log2FoldChange, color = Phylum)) + 
     geom_point(size=5) + 
     scale_color_manual(values=c("#45778f","#802f31"))+ 
     theme(axis.text.x = element_text(angle = -70,hjust=0.06, vjust=.8,size=14, 
             color = "black",face = "italic"),
        axis.text.y = element_text(size=16, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=.7), 
        text = element_text(size=18));PhylumSpp



#GRAPH 3------------------------------------------------------------------------------------------

#Genus to species ...................................................................

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




# PLOT GENUS TO SPECIES GRAPH ........................................................
colors<-c("#69999e","#9888a5","#f7d2f7","#f4df48","#1092b6","#22b14d",
          "#4f3773","#e9d7b4","#38381c","#6cb5f4","#e6ba0b","#8a3500",
          "#2755b1","#676c33", "#5b1700","#f1d68d","#a4b7a2","#b0bac4",
          "#445b39","#e7661b","#c48211")
          
        
GenSpp<-ggplot(sigtabgen, aes(x = Species, y = log2FoldChange, color = Genus)) +
            geom_point(size=5) + 
            scale_color_manual(values=colors)+ 
            theme(axis.text.x = element_text(angle = -70,hjust=0.06, vjust=.8,size=12),
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
pdf("Analysis/Deseq/Ectomycorrhizal/Graphs/Phylum-Genus.pdf", height=6, width=8)
PhylumGen
dev.off()

pdf("Analysis/Deseq/Ectomycorrhizal/Graphs/Phylum-Spp.pdf", height=6, width=8)
PhylumSpp
dev.off()

pdf("Analysis/Deseq/Ectomycorrhizal/Graphs/Genus-Spp.pdf", height=8, width=10)
GenSpecies
dev.off()



#
##
###
#************************************************************************************************----
# CREATE AND EXPORT PANELS --------------------------------------------------------------------------
#************************************************************************************************----
#Ran the codes above individulaly for EcM fungi and Saprobes using their respective codes 
#while does results were loaded on R environment, we ran the codes below to create the panels 
#of both EcM and Saprobes for publication


library(patchwork)


DeseqBothGen<-(PhylumGen/PhylumGenSap);DeseqBothGen 
DeseqBothSpp<- (PhylumSpp/PhylumSppSap);DeseqBothSpp


#EXPORT GRAPHS....... ..........................................................................
pdf("Analysis/Deseq/Publication/Panels-EcM-Sap-PhylumGenus.pdf", height=11, width=11)
DeseqBothGen
dev.off()

pdf("Analysis/Deseq/Publication/Panels-EcM-Sap-PhylumSpecies.pdf", height=11, width=9)
  DeseqBothSpp
dev.off()





