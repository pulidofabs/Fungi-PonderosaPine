#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO/")


#Load librarires----------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
#library(readr)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)
library(RAM)

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
physeqBurnSap<-subset_samples(physeqSap,Treatment=="Burned")
sample_names(physeqBurnSap);physeqBurnSap #106 total

physeqUnburnSap<-subset_samples(physeqSap,Treatment=="Unburned")
sample_names(physeqUnburnSap);physeqUnburnSap #106 total


#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----
#*********************************************************************************************************-----

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


#SET THE REFERENCE LEVEL FOR ANALYSIS(COMPARISON)................................................
attach(RelGenTrt1)
RelGenTrt1$Treatment<-as.factor(RelGenTrt1$Treatment);levels(RelGenTrt1$Treatment)
RelGenTrt1$Treatment <- try(relevel(RelGenTrt1$Treatment , "Unburned"));levels(RelGenTrt1$Treatment)
levels(RelGenTrt1$TimeSinceFire)#looks right
detach(RelGenTrt1)


#************************************************************************************************----
#SUBSET DATA ------------------------------------------------------------------------------------
##***********************************************************************************************----
Archaeorhizomyces<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Archaeorhizomyces"), ]
Basidioascus<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Basidioascus"), ]
Calyptrozyma<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Calyptrozyma"), ]
Cladophialophora<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Cladophialophora"), ]
Geminibasidium<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Geminibasidium"), ]
Hygrocybe<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Hygrocybe"), ]
Mortierella<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Mortierella"), ]
Penicillium<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Penicillium"), ]
Pseudeurotium<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Pseudeurotium"), ]
Rasamsonia<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Rasamsonia"), ]
Sagenomella<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Sagenomella"), ]
Solicoccozyma<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Solicoccozyma"), ]
Umbelopsis<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "Umbelopsis"), ]
#Other<-RelGenTrt1sap[which(RelGenTrt1sap$Genus == "< 2% abund."), ]


#*************************************************************************************************************************************************----
#calculate significance in change per treatment ------------------------------------------------------------------------------------------------------
#*************************************************************************************************************************************************----

#Treatment significance per taxa -----------------------------------------------------
Archaeorhizomyces.Kw<-kruskal.test(Archaeorhizomyces$Abundance ~ Archaeorhizomyces$Treatment, data=Archaeorhizomyces);Archaeorhizomyces.Kw #0.0002525

Basidioascus.Kw<-kruskal.test(Basidioascus$Abundance~Basidioascus$Treatment, data=Basidioascus);Basidioascus.Kw # 2.2e-16

Calyptrozyma.Kw<-kruskal.test(Calyptrozyma$Abundance~Calyptrozyma$Treatment, data=Calyptrozyma);Calyptrozyma.Kw # 2.2e-16

Cladophialophora.Kw<-kruskal.test(Cladophialophora$Abundance~Cladophialophora$Treatment, data=Cladophialophora);Cladophialophora.Kw # 2.2e-167

Geminibasidium.Kw<-kruskal.test(Geminibasidium$Abundance~Geminibasidium$Treatment, data=Geminibasidium);Geminibasidium.Kw #0.2788

Hygrocybe.Kw<-kruskal.test(Hygrocybe$Abundance~Hygrocybe$Treatment, data=Hygrocybe);Hygrocybe.Kw # 2.998e-05

Mortierella.Kw<-kruskal.test(Mortierella$Abundance ~ Mortierella$Treatment, data=Mortierella);Mortierella.Kw #9.577e-09

Penicillium.Kw<-kruskal.test(Penicillium$Abundance~Penicillium$Treatment, data=Penicillium);Penicillium.Kw#0.1504

Pseudeurotium.Kw<-kruskal.test(Pseudeurotium$Abundance~Pseudeurotium$Treatment, data=Pseudeurotium);Pseudeurotium.Kw #1.174e-10

Rasamsonia.Kw<-kruskal.test(Rasamsonia$Abundance~Rasamsonia$Treatment, data=Rasamsonia);Rasamsonia.Kw #8.216e-12

Sagenomella.Kw<-kruskal.test(Sagenomella$Abundance~Sagenomella$Treatment, data=Sagenomella);Sagenomella.Kw #0.1268

Solicoccozyma.Kw<-kruskal.test(Solicoccozyma$Abundance~Solicoccozyma$Treatment, data=Solicoccozyma);Solicoccozyma.Kw #0.002757

Umbelopsis.KW<-kruskal.test(Umbelopsis$Abundance~Umbelopsis$Treatment, data=Umbelopsis);Umbelopsis.KW #0.3161




#*********************************************************************************************************************************----
# CREATE A TABLE FOR EASY EXPORT -----------------------------------------------------------------------------------------------------
#*********************************************************************************************************************************----

KruskalWallisTrt<-cbind(Archaeorhizomyces.Kw,Basidioascus.Kw,Calyptrozyma.Kw,Cladophialophora.Kw,
                        Geminibasidium.Kw,Hygrocybe.Kw,Mortierella.Kw,Penicillium.Kw,Pseudeurotium.Kw,
                        Rasamsonia.Kw,Sagenomella.Kw,Solicoccozyma.Kw,Umbelopsis.KW);KruskalWallisTrt

#Export table

dir.create("Analysis/RelativeAbundance/Saprobes/Tables/Treatment")
capture.output(KruskalWallisTrt, file="Analysis/RelativeAbundance/Saprobes/Tables/Treatment/KruskalWallis-Taxa-Trt.csv")





#
#***********************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE -BURNED SAMPLES--------------------------------------------------------------------
##**********************************************************************************************************************----
#
attach(RelGenTsf1Bsap)
AnthracobiaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Anthracobia"), ]
ArchaeorhizomycesB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Archaeorhizomyces"), ]
BasidioascusB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Basidioascus"), ]
CalyptrozymaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Calyptrozyma"), ]
ClavariaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Clavaria"), ]
ConiochaetaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Coniochaeta"), ]
GeminibasidiumB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Geminibasidium"), ]
HoltermanniellaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Holtermanniella"), ]
HygrocybeB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Hygrocybe"), ]
LeohumicolaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Leohumicola"), ]
LeucosporidiumB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Leucosporidium"), ]
MortierellaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Mortierella"), ]
NeurosporaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Neurospora"), ]
OchrocladosporiumB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Ochrocladosporium"), ]
PenicilliumB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Penicillium"), ]
PholiotaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Pholiota"), ]
PseudeurotiumB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Pseudeurotium"), ]
PyronemaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Pyronema"), ]
RasamsoniaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Rasamsonia"), ]
RutstroemiaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Rutstroemia"), ]
SagenomellaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Sagenomella"), ]
SolicoccozymaB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Solicoccozyma"), ]
TephrocybeB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Tephrocybe"), ]
UmbelopsisB<-RelGenTsf1Bsap[which(RelGenTsf1Bsap$Genus == "Umbelopsis"), ]
detach(RelGenTsf1Bsap)


#********************************************************************************************************************************-----
#- SIGNIFICANCE PER TIME SINCE FIRE -------------------------------------------------------------------------------------------------
#********************************************************************************************************************************-----

AnthracobiaGLM<-glm(AnthracobiaB$Abundance ~ AnthracobiaB$TimeSinceFire + (1/Plot), data = AnthracobiaB)
AnthracobiaAOV<-aov(AnthracobiaGLM);summary(AnthracobiaAOV) #1.15e-05 ***
AnthracobiaTK<-TukeyHSD(AnthracobiaAOV); AnthracobiaTK #All in comparison to 3

ArchaeorhizomycesGLM<-glm(ArchaeorhizomycesB$Abundance ~ ArchaeorhizomycesB$TimeSinceFire + (1/Plot), data = ArchaeorhizomycesB)
ArchaeorhizomycesAOV<-aov(ArchaeorhizomycesGLM);summary(ArchaeorhizomycesAOV) #1.28e-11 ***
ArchaeorhizomycesTK<-TukeyHSD(ArchaeorhizomycesAOV); ArchaeorhizomycesTK #All in comparison to 11

BasidioascusGLM<-glm(BasidioascusB$Abundance ~ BasidioascusB$TimeSinceFire + (1/Plot), data = BasidioascusB)
BasidioascusAOV<-aov(BasidioascusGLM);summary(BasidioascusAOV) #0.00866**
BasidioascusTK<-TukeyHSD(BasidioascusAOV); BasidioascusTK# All compared to 11

CalyptrozymaGLM<-glm(CalyptrozymaB$Abundance ~ CalyptrozymaB$TimeSinceFire + (1/Plot), data = CalyptrozymaB)
CalyptrozymaAOV<-aov(CalyptrozymaGLM);summary(CalyptrozymaAOV) #0.0127 *
CalyptrozymaTK<-TukeyHSD(CalyptrozymaAOV); CalyptrozymaTK #only 5-3

ClavariaGLM<-glm(ClavariaB$Abundance ~ ClavariaB$TimeSinceFire + (1/Plot), data = ClavariaB)
ClavariaAOV<-aov(ClavariaGLM);summary(ClavariaAOV) #3.95e-07 ***
ClavariaTK<-TukeyHSD(ClavariaAOV); ClavariaTK#-- all sign except 11-3

ConiochaetaGLM<-glm(ConiochaetaB$Abundance ~ ConiochaetaB$TimeSinceFire + (1/Plot), data = ConiochaetaB)
ConiochaetaAOV<-aov(ConiochaetaGLM);summary(ConiochaetaAOV) #0.0369
ConiochaetaTK<-TukeyHSD(ConiochaetaAOV); ConiochaetaTK #None

GeminibasidiumGLM<-glm(GeminibasidiumB$Abundance ~ GeminibasidiumB$TimeSinceFire + (1/Plot), data = GeminibasidiumB)
GeminibasidiumAOV<-aov(GeminibasidiumGLM);summary(GeminibasidiumAOV) #0.000128 ***
GeminibasidiumTK<-TukeyHSD(GeminibasidiumAOV); GeminibasidiumTK#11-3;5-3;11-2;5-2

HoltermanniellaGLM<-glm(HoltermanniellaB$Abundance ~ HoltermanniellaB$TimeSinceFire + (1/Plot), data = HoltermanniellaB)
HoltermanniellaAOV<-aov(HoltermanniellaGLM);summary(HoltermanniellaAOV) #2.83e-05 ***
HoltermanniellaTK<-TukeyHSD(HoltermanniellaAOV); HoltermanniellaTK#11-5;5-3;5-2

HygrocybeGLM<-glm(HygrocybeB$Abundance ~ HygrocybeB$TimeSinceFire + (1/Plot), data = HygrocybeB)
HygrocybeAOV<-aov(HygrocybeGLM);summary(HygrocybeAOV) #0.124
HygrocybeTK<-TukeyHSD(HygrocybeAOV); HygrocybeTK 

LeohumicolaGLM<-glm(LeohumicolaB$Abundance ~ LeohumicolaB$TimeSinceFire + (1/Plot), data = LeohumicolaB)
LeohumicolaAOV<-aov(LeohumicolaGLM);summary(LeohumicolaAOV) #0.0158 *
LeohumicolaTK<-TukeyHSD(LeohumicolaAOV); LeohumicolaTK #only 11-3

LeucosporidiumGLM<-glm(LeucosporidiumB$Abundance ~ LeucosporidiumB$TimeSinceFire + (1/Plot), data = LeucosporidiumB)
LeucosporidiumAOV<-aov(LeucosporidiumGLM);summary(LeucosporidiumAOV) #1.03e-05 ***
LeucosporidiumTK<-TukeyHSD(LeucosporidiumAOV); LeucosporidiumTK #All to 11

MortierellaGLM<-glm(MortierellaB$Abundance ~ MortierellaB$TimeSinceFire + (1/Plot), data = MortierellaB)
MortierellaAOV<-aov(MortierellaGLM);summary(MortierellaAOV)#0.00113 **
MortierellaTK<-TukeyHSD(MortierellaAOV); MortierellaTK #3-2

NeurosporaGLM<-glm(NeurosporaB$Abundance ~ NeurosporaB$TimeSinceFire + (1/Plot), data = NeurosporaB)
NeurosporaAOV<-aov(NeurosporaGLM);summary(NeurosporaAOV)#0.0606.
NeurosporaTK<-TukeyHSD(NeurosporaAOV); NeurosporaTK 

OchrocladosporiumGLM<-glm(OchrocladosporiumB$Abundance ~ OchrocladosporiumB$TimeSinceFire + (1/Plot), data = OchrocladosporiumB)
OchrocladosporiumAOV<-aov(OchrocladosporiumGLM);summary(OchrocladosporiumAOV)#1.21e-07 ***
OchrocladosporiumTK<-TukeyHSD(OchrocladosporiumAOV); OchrocladosporiumTK#All in compariso to 2

PenicilliumGLM<-glm(PenicilliumB$Abundance ~ PenicilliumB$TimeSinceFire + (1/Plot), data = PenicilliumB)
PenicilliumAOV<-aov(PenicilliumGLM);summary(PenicilliumAOV)#1.55e-06 ***
PenicilliumTK<-TukeyHSD(PenicilliumAOV); PenicilliumTK #All in comparison to 3

PholiotaGLM<-glm(PholiotaB$Abundance ~ PholiotaB$TimeSinceFire + (1/Plot), data = PholiotaB)
PholiotaAOV<-aov(PholiotaGLM);summary(PholiotaAOV)#0.0311 *
PholiotaTK<-TukeyHSD(PholiotaAOV); PholiotaTK#

PseudeurotiumGLM<-glm(PseudeurotiumB$Abundance ~ PseudeurotiumB$TimeSinceFire + (1/Plot), data = PseudeurotiumB)
PseudeurotiumAOV<-aov(PseudeurotiumGLM);summary(PseudeurotiumAOV) #0.000951 ***
PseudeurotiumTK<-TukeyHSD(PseudeurotiumAOV); PseudeurotiumTK# All years in comparison to 3

PyronemaGLM<-glm(PyronemaB$Abundance ~ PyronemaB$TimeSinceFire + (1/Plot), data = PyronemaB)
PyronemaAOV<-aov(PyronemaGLM);summary(PyronemaAOV) #0.00142**
PyronemaTK<-TukeyHSD(PyronemaAOV); PyronemaTK# 11-2

RasamsoniaGLM<-glm(RasamsoniaB$Abundance ~ RasamsoniaB$TimeSinceFire + (1/Plot), data = RasamsoniaB)
RasamsoniaAOV<-aov(RasamsoniaGLM);summary(RasamsoniaAOV) #0.000787 ***
RasamsoniaTK<-TukeyHSD(RasamsoniaAOV); RasamsoniaTK#all except 11-3; 5-3

RutstroemiaGLM<-glm(RutstroemiaB$Abundance ~ RutstroemiaB$TimeSinceFire + (1/Plot), data = RutstroemiaB)
RutstroemiaAOV<-aov(RutstroemiaGLM);summary(RutstroemiaAOV) #0.485
RutstroemiaTK<-TukeyHSD(RutstroemiaAOV); RutstroemiaTK#

SagenomellaGLM<-glm(SagenomellaB$Abundance ~ SagenomellaB$TimeSinceFire + (1/Plot), data = SagenomellaB)
SagenomellaAOV<-aov(SagenomellaGLM);summary(SagenomellaAOV) #2.49e-07 ***
SagenomellaTK<-TukeyHSD(SagenomellaAOV); SagenomellaTK#All except 11-5; 3-2

SolicoccozymaGLM<-glm(SolicoccozymaB$Abundance ~ SolicoccozymaB$TimeSinceFire + (1/Plot), data = SolicoccozymaB)
SolicoccozymaAOV<-aov(SolicoccozymaGLM);summary(SolicoccozymaAOV) #0.0118*
SolicoccozymaTK<-TukeyHSD(SolicoccozymaAOV); SolicoccozymaTK #5-2;11-5

TephrocybeGLM<-glm(TephrocybeB$Abundance ~ TephrocybeB$TimeSinceFire + (1/Plot), data = TephrocybeB)
TephrocybeAOV<-aov(TephrocybeGLM);summary(TephrocybeAOV) #0.0119*
TephrocybeTK<-TukeyHSD(TephrocybeAOV); TephrocybeTK #none

UmbelopsisGLM<-glm(UmbelopsisB$Abundance ~ UmbelopsisB$TimeSinceFire + (1/Plot), data = UmbelopsisB)
UmbelopsisAOV<-aov(UmbelopsisGLM);summary(UmbelopsisAOV) #2.06e-09 ***
UmbelopsisTK<-TukeyHSD(UmbelopsisAOV); UmbelopsisTK #All in comparison except 3-2; 11-5



#Create a table and export the results------------------------------------------------------------------------------------------------
dir.create("Analysis/RelativeAbundance/Saprobes/Tables/TSFburn")

capture.output(summary(AnthracobiaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Anthracobia-glm.csv")
capture.output(summary(ArchaeorhizomycesAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Archaeorhizomyces-glm.csv")
capture.output(summary(BasidioascusAOV),file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Basidioascus-glm.csv")
capture.output(summary(CalyptrozymaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Calyptrozyma-glm.csv")
capture.output(summary(ClavariaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Clavaria-glm.csv")
capture.output(summary(ConiochaetaAOV),file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Coniochaeta-glm.csv")
capture.output(summary(GeminibasidiumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Geminibasidium-glm.csv")
capture.output(summary(HoltermanniellaAOV), file = "Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Holtermanniella-glm.csv")
capture.output(summary(HygrocybeAOV),file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Hygrocybe-glm.csv")
capture.output(summary(LeohumicolaAOV), file = "Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Leohumicola-glm.csv")
capture.output(summary(LeucosporidiumAOV),file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Leucosporidium-glm.csv")
capture.output(summary(MortierellaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Mortierella-glm.csv")
capture.output(summary(NeurosporaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Neurospora-glm.csv")
capture.output(summary(OchrocladosporiumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Ochrocladosporium-glm.csv")
capture.output(summary(PenicilliumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Penicillium-glm.csv")
capture.output(summary(PholiotaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Pholiota-glm.csv")
capture.output(summary(PseudeurotiumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Pseudeurotium-glm.csv")
capture.output(summary(PyronemaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Pyronema-glm.csv")
capture.output(summary(RasamsoniaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Rasamsonia-glm.csv")
capture.output(summary(RutstroemiaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Rutstroemia-glm.csv")
capture.output(summary(SagenomellaAOV),file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Sagenomella-glm.csv")
capture.output(summary(SolicoccozymaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Solicoccozyma-glm.csv")
capture.output(summary(TephrocybeAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Tephrocybe-glm.csv")
capture.output(summary(UmbelopsisAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/Umbelopsis-glm.csv")



#Export tukeys test
capture.output(AnthracobiaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/AnthracobiaTK.csv")
capture.output(ArchaeorhizomycesTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/ArchaeorhizomycesTK.csv")
capture.output(BasidioascusTK,file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/BasidioascusTK.csv")
capture.output(CalyptrozymaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/CalyptrozymaTK.csv")
capture.output(ClavariaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/ClavariaTK.csv")
capture.output(ConiochaetaTK,file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/ConiochaetaTK.csv")
capture.output(GeminibasidiumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/GeminibasidiumTK.csv")
capture.output(HoltermanniellaTK, file = "Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/HoltermanniellaTK.csv")
capture.output(HygrocybeTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/HygrocybeTK.csv")
capture.output(LeohumicolaTK, file = "Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/LeohumicolaTK.csv")
capture.output(LeucosporidiumTK,file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/LeucosporidiumTK.csv")
capture.output(MortierellaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/MortierellaTK.csv")
capture.output(NeurosporaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/NeurosporaTK.csv")
capture.output(OchrocladosporiumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/OchrocladosporiumTK.csv")
capture.output(PenicilliumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/PenicilliumTK.csv")
capture.output(PholiotaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/PholiotaTK.csv")
capture.output(PseudeurotiumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/PseudeurotiumTK.csv")
capture.output(PyronemaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/PyronemaTK.csv")
capture.output(RasamsoniaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/RasamsoniaTK.csv")
capture.output(RutstroemiaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/RutstroemiaTK.csv")
capture.output(SagenomellaTK,file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/SagenomellaTK.csv")
capture.output(SolicoccozymaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/SolicoccozymaTK.csv")
capture.output(TephrocybeTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/TephrocybeTK.csv")
capture.output(UmbelopsisTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFburn/UmbelopsisTK.csv")



#*************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE - UNBURNED--------------------------------------------------------------------------
##************************************************************************************************************************----
#
attach(RelGenTsf1UnSap)

ArchaeorhizomycesUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Archaeorhizomyces"), ]
CladophialophoraUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Cladophialophora"), ]
ClavariaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Clavaria"), ]
ConiochaetaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Coniochaeta"), ]
GeminibasidiumUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Geminibasidium"), ]
HumicolaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Humicola"), ]
HygrocybeUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Hygrocybe"), ]
MortierellaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Mortierella"), ]
PenicilliumUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Penicillium"), ]
PseudeurotiumUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Pseudeurotium"), ]
SagenomellaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Sagenomella"), ]
SolicoccozymaUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Solicoccozyma"), ]
UmbelopsisUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "Umbelopsis"), ]
unidentifiedUn<-RelGenTsf1UnSap[which(RelGenTsf1UnSap$Genus == "unidentified"), ]





#*****************************************************************************************************************************************-----
#- SIGNIFICANCE PER TIME SINCE FIRE -----------------------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************-----
ArchaeorhizomycesGLM<-glm(ArchaeorhizomycesUn$Abundance ~ ArchaeorhizomycesUn$TimeSinceFire + (1/Plot), data = ArchaeorhizomycesUn)
ArchaeorhizomycesAOV<-aov(ArchaeorhizomycesGLM);summary(ArchaeorhizomycesAOV) #0.000668 ***
ArchaeorhizomycesTK<-TukeyHSD(ArchaeorhizomycesAOV); ArchaeorhizomycesTK #only 11-3;5-3;5-2

CladophialophoraGLM<-glm(CladophialophoraUn$Abundance ~ CladophialophoraUn$TimeSinceFire + (1/Plot), data = CladophialophoraUn)
CladophialophoraAOV<-aov(CladophialophoraGLM);summary(CladophialophoraAOV) #0.078.
CladophialophoraTK<-TukeyHSD(CladophialophoraAOV); CladophialophoraTK

ClavariaGLM<-glm(ClavariaUn$Abundance ~ ClavariaUn$TimeSinceFire + (1/Plot), data = ClavariaUn)
ClavariaAOV<-aov(ClavariaGLM);summary(ClavariaAOV) #0.00579 **
ClavariaTK<-TukeyHSD(ClavariaAOV); ClavariaTK #All to 3

ConiochaetaGLM<-glm(ConiochaetaUn$Abundance ~ ConiochaetaUn$TimeSinceFire + (1/Plot), data = ConiochaetaUn)
ConiochaetaAOV<-aov(ConiochaetaGLM);summary(ConiochaetaAOV) #7.16e-07 ***
ConiochaetaTK<-TukeyHSD(ConiochaetaAOV); ConiochaetaTK #only 5-2;11-2;5-3

GeminibasidiumGLM<-glm(GeminibasidiumUn$Abundance ~ GeminibasidiumUn$TimeSinceFire + (1/Plot), data = GeminibasidiumUn)
GeminibasidiumAOV<-aov(GeminibasidiumGLM);summary(GeminibasidiumAOV) #1.13e-06 ***
GeminibasidiumTK<-TukeyHSD(GeminibasidiumAOV); GeminibasidiumTK #All except 5-3;5-2;3-2

HumicolaGLM<-glm(HumicolaUn$Abundance ~ HumicolaUn$TimeSinceFire + (1/Plot), data = HumicolaUn)
HumicolaAOV<-aov(HumicolaGLM);summary(HumicolaAOV) #3.68e-10 ***
HumicolaTK<-TukeyHSD(HumicolaAOV); HumicolaTK #All in comparison to 5

HygrocybeGLM<-glm(HygrocybeUn$Abundance ~ HygrocybeUn$TimeSinceFire + (1/Plot), data = HygrocybeUn)
HygrocybeAOV<-aov(HygrocybeGLM);summary(HygrocybeAOV) #7.39e-06 ***
HygrocybeTK<-TukeyHSD(HygrocybeAOV); HygrocybeTK #only 5-2;11-2;5-2

MortierellaGLM<-glm(MortierellaUn$Abundance ~ MortierellaUn$TimeSinceFire + (1/Plot), data = MortierellaUn)
MortierellaAOV<-aov(MortierellaGLM);summary(MortierellaAOV) #1.32e-08 ***
MortierellaTK<-TukeyHSD(MortierellaAOV); MortierellaTK #All except 3-2

PenicilliumGLM<-glm(PenicilliumUn$Abundance ~ PenicilliumUn$TimeSinceFire + (1/Plot), data = PenicilliumUn)
PenicilliumAOV<-aov(PenicilliumGLM);summary(PenicilliumAOV) #7.59e-12 ***
PenicilliumTK<-TukeyHSD(PenicilliumAOV); PenicilliumTK #only 3-2

PseudeurotiumGLM<-glm(PseudeurotiumUn$Abundance ~ PseudeurotiumUn$TimeSinceFire + (1/Plot), data = PseudeurotiumUn)
PseudeurotiumAOV<-aov(PseudeurotiumGLM);summary(PseudeurotiumAOV) #4.64e-09 ***
PseudeurotiumTK<-TukeyHSD(PseudeurotiumAOV); PseudeurotiumTK #only 3-2

SagenomellaGLM<-glm(SagenomellaUn$Abundance ~ SagenomellaUn$TimeSinceFire + (1/Plot), data = SagenomellaUn)
SagenomellaAOV<-aov(SagenomellaGLM);summary(SagenomellaAOV) #6.67e-14 ***
SagenomellaTK<-TukeyHSD(SagenomellaAOV); SagenomellaTK #All except 5-2

SolicoccozymaGLM<-glm(SolicoccozymaUn$Abundance ~ SolicoccozymaUn$TimeSinceFire + (1/Plot), data = SolicoccozymaUn)
SolicoccozymaAOV<-aov(SolicoccozymaGLM);summary(SolicoccozymaAOV) #<2e-16 ***
SolicoccozymaTK<-TukeyHSD(SolicoccozymaAOV); SolicoccozymaTK #All in comparison to 5

UmbelopsisGLM<-glm(UmbelopsisUn$Abundance ~ UmbelopsisUn$TimeSinceFire + (1/Plot), data = UmbelopsisUn)
UmbelopsisAOV<-aov(UmbelopsisGLM);summary(UmbelopsisAOV) #0.00262 **
UmbelopsisTK<-TukeyHSD(UmbelopsisAOV); UmbelopsisTK 

unidentifiedGLM<-glm(unidentifiedUn$Abundance ~ unidentifiedUn$TimeSinceFire + (1/Plot), data = unidentifiedUn)
unidentifiedAOV<-aov(unidentifiedGLM);summary(unidentifiedAOV) #0.00183 **
unidentifiedTK<-TukeyHSD(unidentifiedAOV); unidentifiedTK #only 5-2;11-3;11-5




#****************************************************************************************************************************************----
# ------------------- EXPORT THE RESULTS ----------------------------------------------------------------------------------------------------
#****************************************************************************************************************************************----

#ANOVA SIGNIFICANCE......................................................................................................................
dir.create("Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn")
capture.output(summary(ArchaeorhizomycesAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Archaeorhizomyces-glm.csv")
capture.output(summary(CladophialophoraAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Cladophialophora-glm.csv")
capture.output(summary(ClavariaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Clavaria-glm.csv")
capture.output(summary(ConiochaetaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Coniochaeta-glm.csv")
capture.output(summary(GeminibasidiumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Geminibasidium-glm.csv")
capture.output(summary(HumicolaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Humicola-glm.csv")
capture.output(summary(HygrocybeAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Hygrocybe-glm.csv")
capture.output(summary(MortierellaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Mortierella-glm.csv")
capture.output(summary(PenicilliumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Penicillium-glm.csv")
capture.output(summary(PseudeurotiumAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Pseudeutorium-glm.csv")
capture.output(summary(SagenomellaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Sagenomella-glm.csv")
capture.output(summary(SolicoccozymaAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Solicoccozyma-glm.csv")
capture.output(summary(UmbelopsisAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Umbelopsis-glm.csv")
capture.output(summary(unidentifiedAOV), file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/unidentified-glm.csv")


#TUKEY SIGNIFICANCE......................................................................................................................
capture.output(ArchaeorhizomycesTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Archaeorhizomyces-Tukey.csv")
capture.output(CladophialophoraTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Cladophialophora-Tukey.csv")
capture.output(ClavariaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Clavaria-Tukey.csv")
capture.output(ConiochaetaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Coniochaeta-Tukey.csv")
capture.output(GeminibasidiumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Geminibasidium-Tukey.csv")
capture.output(HumicolaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Humicola-Tukey.csv")
capture.output(HygrocybeTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Hygrocybe-Tukey.csv")
capture.output(MortierellaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Mortierella-Tukey.csv")
capture.output(PenicilliumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Penicillium-Tukey.csv")
capture.output(PseudeurotiumTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Pseudeutorium-Tukey.csv")
capture.output(SagenomellaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Sagenomella-Tukey.csv")
capture.output(SolicoccozymaTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Solicoccozyma-Tukey.csv")
capture.output(UmbelopsisTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/Umbelopsis-Tukey.csv")
capture.output(unidentifiedTK, file="Analysis/RelativeAbundance/Saprobes/Tables/TSFunburn/unidentified-Tukey.csv")




















]
