
#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/PIPO/")


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


#***************************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#***************************************************************************************************************************----
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


#***********************************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************************-----

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


#Relevel data to appropiate reference level ...........................................................
attach(RelGenTrt1)
RelGenTrt1$Treatment<-as.factor(RelGenTrt1$Treatment);levels(RelGenTrt1$Treatment)
RelGenTrt1$Treatment <- try(relevel(RelGenTrt1$Treatment , "Unburned"));levels(RelGenTrt1$Treatment)
levels(RelGenTrt1$TimeSinceFire)#looks right
detach(RelGenTrt1)


#EXPORT FILES .............................................................................
write.csv(RelGenTrt1, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundaceTrt-New.csv")
write.csv(RelGenTsf1B, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundace-TSF-Burn-New.csv")
write.csv(RelGenTsf1Un, "Analysis/RelativeAbundance/Ectomycorrhizal/Tables/RelativeAbundace-TSF-Unburn-New.csv")


#***********************************************************************************************************************-----
# SUBSET DATA PER TAXA OF INTEREST -------------------------------------------------------------------------------
#***********************************************************************************************************************-----


#SUBSET THE DATA---------------------------------------------------------------------
Cenococcum<-RelGenTrt1[which(RelGenTrt1$Genus == "Cenococcum"), ]
Chromelosporium<-RelGenTrt1[which(RelGenTrt1$Genus == "Chromelosporium"), ]
Cortinarius<-RelGenTrt1[which(RelGenTrt1$Genus == "Cortinarius"), ]
Helvellosebacina<-RelGenTrt1[which(RelGenTrt1$Genus == "Helvellosebacina"), ]
Hygrophorus<-RelGenTrt1[which(RelGenTrt1$Genus == "Hygrophorus"), ]
Inocybe<-RelGenTrt1[which(RelGenTrt1$Genus == "Inocybe"), ]
Lyophyllum<-RelGenTrt1[which(RelGenTrt1$Genus == "Lyophyllum"), ]
Peziza<-RelGenTrt1[which(RelGenTrt1$Genus == "Peziza"), ]
Piloderma<-RelGenTrt1[which(RelGenTrt1$Genus == "Piloderma"), ]
Pustularia<-RelGenTrt1[which(RelGenTrt1$Genus == "Pustularia"), ]
Rhizopogon<-RelGenTrt1[which(RelGenTrt1$Genus == "Rhizopogon"), ]
Russula<-RelGenTrt1[which(RelGenTrt1$Genus == "Russula"), ]
Thelephora<-RelGenTrt1[which(RelGenTrt1$Genus == "Thelephora"), ]
Tomentella<-RelGenTrt1[which(RelGenTrt1$Genus == "Tomentella"), ]
Tuber<-RelGenTrt1[which(RelGenTrt1$Genus == "Tuber"), ]
Wilcoxina<-RelGenTrt1[which(RelGenTrt1$Genus == "Wilcoxina"), ]
#Other<-RelGenTrt1[which(RelGenTrt1$Genus == "< 2% abund."), ]




#*****************************************************************************************************************************************----
#SIGNIFICANCE PER TAXA PER TREATMENT----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----

#Treatment significance per taxa -----------------------------------------------------
Cenococcum.Kw<-kruskal.test(Cenococcum$Abundance ~ Cenococcum$Treatment, data=Cenococcum);Cenococcum.Kw #0.2999

Chromelosporium.Kw<-kruskal.test(Chromelosporium$Abundance ~ Chromelosporium$Treatment, data=Chromelosporium);Chromelosporium.Kw #2.2e-16

Cortinarius.Kw<-kruskal.test(Cortinarius$Abundance~Cortinarius$Treatment, data=Cortinarius);Cortinarius.Kw #0.001288

Helvellosebacina.Kw<-kruskal.test(Helvellosebacina$Abundance~Helvellosebacina$Treatment, data=Helvellosebacina);Helvellosebacina.Kw #0.001741

Hygrophorus.Kw<-kruskal.test(Hygrophorus$Abundance~Hygrophorus$Treatment, data=Hygrophorus);Hygrophorus.Kw #1.089e-05

Inocybe.Kw<-kruskal.test(Inocybe$Abundance~Inocybe$Treatment, data=Inocybe);Inocybe.Kw #1.101e-07

Lyophyllum.Kw<-kruskal.test(Lyophyllum$Abundance~Lyophyllum$Treatment, data=Lyophyllum);Lyophyllum.Kw #< 2.2e-16

Peziza.Kw<-kruskal.test(Peziza$Abundance ~ Peziza$Treatment, data=Peziza);Peziza.Kw #1.01e-05

Piloderma.Kw<-kruskal.test(Piloderma$Abundance~Piloderma$Treatment, data=Piloderma);Piloderma.Kw#9.223e-07

Pustularia.Kw<-kruskal.test(Pustularia$Abundance~Pustularia$Treatment, data=Pustularia);Pustularia.Kw #1.73e-11

Rhizopogon.Kw<-kruskal.test(Rhizopogon$Abundance~Rhizopogon$Treatment, data=Rhizopogon);Rhizopogon.Kw #0.0001147

Russula.Kw<-kruskal.test(Russula$Abundance~Russula$Treatment, data=Russula);Russula.Kw #< 2.2e-16

Thelephora.Kw<-kruskal.test(Thelephora$Abundance~Thelephora$Treatment, data=Thelephora);Thelephora.Kw #4.047e-06

Tomentella.KW<-kruskal.test(Tomentella$Abundance~Tomentella$Treatment, data=Tomentella);Tomentella.KW #0.0001619

Tuber.KW<-kruskal.test(Tuber$Abundance~Tuber$Treatment, data=Tuber);Tuber.KW #0.01994

Wilcox.KW<-kruskal.test(Wilcoxina$Abundance~Wilcoxina$Treatment, data=Wilcoxina);Wilcox.KW #0.1596



#BIND RESULTS FOR EASY EXPORT ----------------------------------------------------------------------------------------
KruskalWallisTrt<-cbind(Cenococcum.Kw, Chromelosporium.Kw,Cortinarius.Kw,
                  Helvellosebacina.Kw,Hygrophorus.Kw,Inocybe.Kw,Lyophyllum.Kw,
                  Peziza.Kw,Piloderma.Kw,Pustularia.Kw,Rhizopogon.Kw,Russula.Kw,
                  Thelephora.Kw,Tomentella.KW,Tuber.KW,Wilcox.KW);KruskalWallisTrt



#Export table.....................................................................................................
capture.output(KruskalWallisTrt, 
               file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/Treatment/KruskalWallis-Taxa-Trt.csv")




##
##
#************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE -BURNED Ectomycorrhizal ----------------------------------------------------------
##***********************************************************************************************************************----
#
attach(RelGenTsf1B)

#SUBSET DATA BY TAXA OF INTEREST ....................................................
CenococcumB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Cenococcum"), ]
ChromelosporiumB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Chromelosporium"), ]
CortinariusB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Cortinarius"), ]
GeopyxisB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Geopyxis"), ]
HebelomaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Hebeloma"), ]
InocybeB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Inocybe"), ]
LyophyllumB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Lyophyllum"), ]
MorchellaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Morchella"), ]
PezizaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Peziza"), ]
PustulariaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Pustularia"), ]
RhizopogonB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Rhizopogon"), ]
SistotremaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Sistotrema"), ]
ThelephoraB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Thelephora"), ]
TrichophaeaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Trichophaea"), ]
WilcoxinaB<-RelGenTsf1B[which(RelGenTsf1B$Genus == "Wilcoxina"), ]



##
##
#****************************************************************************************************************************-----
#- SIGNIFICANCE PER TIME SINCE FIRE FOR BURNED DATA ------------------------------------------------------------------------------
#****************************************************************************************************************************-----

CenococcumGLM<-glm(CenococcumB$Abundance ~ CenococcumB$TimeSinceFire + (1/Plot), data = CenococcumB)
CenoAOV<-aov(CenococcumGLM);summary(CenoAOV)#0.0334 *
CenoTK<-TukeyHSD(CenoAOV); CenoTK #only 5-1; 11-2

ChromelosporiumGLM<-glm(ChromelosporiumB$Abundance ~ ChromelosporiumB$TimeSinceFire + (1/Plot), data = ChromelosporiumB)
ChromAOV<-aov(ChromelosporiumGLM);summary(ChromAOV) #6.17e-07 ***
ChromTK<-TukeyHSD(ChromAOV); ChromTK # 5-2, 5-3, 11-5

CortinariusGLM<-glm(CortinariusB$Abundance ~ CortinariusB$TimeSinceFire + (1/Plot), data = CortinariusB)
CortAOV<-aov(CortinariusGLM);summary(CortAOV) #0.17
CortTK<-TukeyHSD(CortAOV); CortTK

GeopyxisGLM<-glm(GeopyxisB$Abundance ~ GeopyxisB$TimeSinceFire + (1/Plot), data = GeopyxisB)
GeopyxisAOV<-aov(GeopyxisGLM);summary(GeopyxisAOV) #3.22e-06 ***
GeopyxisTK<-TukeyHSD(GeopyxisAOV); GeopyxisTK#-- all sign except 11-3;11-5;5-3

HebelomaGLM<-glm(HebelomaB$Abundance ~ HebelomaB$TimeSinceFire + (1/Plot), data = HebelomaB)
HebelomaAOV<-aov(HebelomaGLM);summary(HebelomaAOV) #0.14
HebelomaTK<-TukeyHSD(HebelomaAOV); HebelomaTK 

InocybeGLM<-glm(InocybeB$Abundance ~ InocybeB$TimeSinceFire + (1/Plot), data = InocybeB)
InocybeAOV<-aov(InocybeGLM);summary(InocybeAOV) #0.0192
InocybeTK<-TukeyHSD(InocybeAOV); InocybeTK#3-2

LyophyllumGLM<-glm(LyophyllumB$Abundance ~ LyophyllumB$TimeSinceFire + (1/Plot), data = LyophyllumB)
LyophyllumAOV<-aov(LyophyllumGLM);summary(LyophyllumAOV) #0.0118
LyophyllumTK<-TukeyHSD(LyophyllumAOV); LyophyllumTK #11-3; 11-5

MorchellaGLM<-glm(MorchellaB$Abundance ~ MorchellaB$TimeSinceFire + (1/Plot), data = MorchellaB)
MorchellaAOV<-aov(MorchellaGLM);summary(MorchellaAOV) #0.253
MorchellaTK<-TukeyHSD(MorchellaAOV); MorchellaTK #

PezizaGLM<-glm(PezizaB$Abundance ~ PezizaB$TimeSinceFire + (1/Plot), data = PezizaB)
PezizaAOV<-aov(PezizaGLM);summary(PezizaAOV)#0.301
PezizaTK<-TukeyHSD(PezizaAOV); PezizaTK

PustulariaGLM<-glm(PustulariaB$Abundance ~ PustulariaB$TimeSinceFire + (1/Plot), data = PustulariaB)
PustulariaAOV<-aov(PustulariaGLM);summary(PustulariaAOV) #4.29e-09 ***
PustulariaTK<-TukeyHSD(PustulariaAOV); PustulariaTK# All years in comparison to 2

RhizopogonGLM<-glm(RhizopogonB$Abundance ~ RhizopogonB$TimeSinceFire + (1/Plot), data = RhizopogonB)
RhizopogonAOV<-aov(RhizopogonGLM);summary(RhizopogonAOV) #0.111
RhizopogonTK<-TukeyHSD(RhizopogonAOV); RhizopogonTK

SistotremaGLM<-glm(SistotremaB$Abundance ~ SistotremaB$TimeSinceFire + (1/Plot), data = SistotremaB)
SistotremaAOV<-aov(SistotremaGLM);summary(SistotremaAOV) #0.524
SistotremaTK<-TukeyHSD(SistotremaAOV); SistotremaTK

ThelephoraGLM<-glm(ThelephoraB$Abundance ~ ThelephoraB$TimeSinceFire + (1/Plot), data = ThelephoraB)
ThelephoraAOV<-aov(ThelephoraGLM);summary(ThelephoraAOV) #0.373
ThelephoraTK<-TukeyHSD(ThelephoraAOV); ThelephoraTK

TrichophaeaGLM<-glm(TrichophaeaB$Abundance ~ TrichophaeaB$TimeSinceFire + (1/Plot), data = TrichophaeaB)
TrichophaeaAOV<-aov(TrichophaeaGLM);summary(TrichophaeaAOV) #0.147
TrichophaeaTK<-TukeyHSD(TrichophaeaAOV); TrichophaeaTK 

WilcoxinaGLM<-glm(WilcoxinaB$Abundance ~ WilcoxinaB$TimeSinceFire + (1/Plot), data = WilcoxinaB)
WilcoxinaAOV<-aov(WilcoxinaGLM);summary(WilcoxinaAOV) #6.29e-07 ***
WilcoxinaTK<-TukeyHSD(WilcoxinaAOV); WilcoxinaTK #All in comparison except 3-2; 5-3



#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create("Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn")

#Export ANOVA values ...............................................................................................................
capture.output(summary(CenoAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Cennococum-glm.csv")
capture.output(summary(ChromAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Chromelosporium-glm.csv")
capture.output(summary(CortAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Cortinarius-glm.csv")
capture.output(summary(GeopyxisAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Geopyxis-glm.csv")
capture.output(summary(HebelomaAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Hebeloma-glm.csv")
capture.output(summary(InocybeAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Inocybe-glm.csv")
capture.output(summary(LyophyllumAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Lyophyllum-glm.csv")
capture.output(summary(MorchellaAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Morchella-glm.csv")
capture.output(summary(PezizaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Peziza-glm.csv")
capture.output(summary(PustulariaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Pustularia-glm.csv")
capture.output(summary(RhizopogonAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Rhizopogon-glm.csv")
capture.output(summary(SistotremaAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Sistotrema-glm.csv")
capture.output(summary(ThelephoraAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Telephora-glm.csv")
capture.output(summary(TrichophaeaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Trichophaea-glm.csv")
capture.output(summary(WilcoxinaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Wilcoxina-glm.csv")



#Export tukeys test...............................................................................................................
capture.output(CenoTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Cennococum-Tukey.csv")
capture.output(ChromTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Chromelosporium-Tukey.csv")
capture.output(CortTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Cortinarius-Tukey.csv")
capture.output(GeopyxisTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Geopyxis-Tukey.csv")
capture.output(HebelomaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Hebeloma-Tukey.csv")
capture.output(InocybeTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Inocybe-Tukey.csv")
capture.output(LyophyllumTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Lyophyllum-Tukey.csv")
capture.output(MorchellaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Morchella-Tukey.csv")
capture.output(PezizaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Peziza-Tukey.csv")
capture.output(PustulariaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Pustularia-Tukey.csv")
capture.output(RhizopogonTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Rhizopogon-Tukey.csv")
capture.output(SistotremaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Sistotrema-Tukey.csv")
capture.output(ThelephoraTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Thelephora-Tukey.csv")
capture.output(TrichophaeaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Trichophaea-Tukey.csv")
capture.output(WilcoxinaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFburn/Wilcoxina-Tukey.csv")

detach(RelGenTsf1B)


#***************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE -UNBURNED Ectomycorrhizal ------------------------------------------------------------
##**************************************************************************************************************************----
#
attach(RelGenTsf1Un)

AmphinemaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Amphinema"), ]
CortinariusUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Cortinarius"), ]
EntolomaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Entoloma"), ]
HelvellosebacinaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Helvellosebacina"), ]
HygrophorusUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Hygrophorus"), ]
InocybeUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Inocybe"), ]
PilodermaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Piloderma"), ]
RhizopogonUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Rhizopogon"), ]
RussulaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Russula"), ]
SistotremaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Sistotrema"), ]
SuillusUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Suillus"), ]
ThelephoraUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Thelephora"), ]
TomentellaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Tomentella"), ]
TricholomaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Tricholoma"), ]
TuberUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Tuber"), ]
TylosporaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Tylospora"), ]
WilcoxinaUn<-RelGenTsf1Un[which(RelGenTsf1Un$Genus == "Wilcoxina"), ]



#***********************************************************************************************************-----
#- SIGNIFICANCE PER TIME SINCE FIRE ----------------------------------------------------------------------------
#***********************************************************************************************************-----
AmphinemaGLM<-glm(AmphinemaUn$Abundance ~ AmphinemaUn$TimeSinceFire + (1/Plot), data = AmphinemaUn)
AmphinemaAOV<-aov(AmphinemaGLM);summary(AmphinemaAOV) #0.116
AmphinemaTK<-TukeyHSD(AmphinemaAOV); AmphinemaTK#

CortinariusGLM<-glm(CortinariusUn$Abundance ~ CortinariusUn$TimeSinceFire + (1/Plot), data = CortinariusUn)
CortAOV<-aov(CortinariusGLM);summary(CortAOV) #0.0159 *
CortTK<-TukeyHSD(CortAOV); CortTK #5-2; 5-3; 

EntolomaGLM<-glm(EntolomaUn$Abundance ~ EntolomaUn$TimeSinceFire + (1/Plot), data = EntolomaUn)
EntolomaAOV<-aov(EntolomaGLM);summary(EntolomaAOV) #0.391
EntolomaTK<-TukeyHSD(EntolomaAOV); EntolomaTK #

HelvellosebacinaGLM<-glm(HelvellosebacinaUn$Abundance ~ HelvellosebacinaUn$TimeSinceFire + (1/Plot), data = HelvellosebacinaUn)
HelvellosebacinaAOV<-aov(HelvellosebacinaGLM);summary(HelvellosebacinaAOV) #0.00177 **
HelvellosebacinaTK<-TukeyHSD(HelvellosebacinaAOV); HelvellosebacinaTK #All years w 5-2;5-3;11-5

HygrophorusGLM<-glm(HygrophorusUn$Abundance ~ HygrophorusUn$TimeSinceFire + (1/Plot), data = HygrophorusUn)
HygrophorusAOV<-aov(HygrophorusGLM);summary(HygrophorusAOV) #0.00265 **
HygrophorusTK<-TukeyHSD(HygrophorusAOV); HygrophorusTK#-- all years compared w 2

InocybeGLM<-glm(InocybeUn$Abundance ~ InocybeUn$TimeSinceFire + (1/Plot), data = InocybeUn)
InocybeAOV<-aov(InocybeGLM);summary(InocybeAOV) #0.00175 **
InocybeTK<-TukeyHSD(InocybeAOV); InocybeTK#Only 5-3

PilodermaGLM<-glm(PilodermaUn$Abundance ~ PilodermaUn$TimeSinceFire + (1/Plot), data = PilodermaUn)
PilodermaAOV<-aov(PilodermaGLM);summary(PilodermaAOV) #7.07e-10 ***
PilodermaTK<-TukeyHSD(PilodermaAOV); PilodermaTK #All compared to 2

RhizopogonGLM<-glm(RhizopogonUn$Abundance ~ RhizopogonUn$TimeSinceFire + (1/Plot), data = RhizopogonUn)
RhizopogonAOV<-aov(RhizopogonGLM);summary(RhizopogonAOV) #0.132
RhizopogonTK<-TukeyHSD(RhizopogonAOV); RhizopogonTK

RussulaGLM<-glm(RussulaUn$Abundance ~ RussulaUn$TimeSinceFire + (1/Plot), data = RussulaUn)
RussulaAOV<-aov(RussulaGLM);summary(RussulaAOV) #1.05e-06 ***
RussulaTK<-TukeyHSD(RussulaAOV); RussulaTK # All in comparison to 5

SistotremaGLM<-glm(SistotremaUn$Abundance ~ SistotremaUn$TimeSinceFire + (1/Plot), data = SistotremaUn)
SistotremaAOV<-aov(SistotremaGLM);summary(SistotremaAOV)#1.31e-05 ***
SistotremaTK<-TukeyHSD(SistotremaAOV); SistotremaTK #All in comparison to 3

SuillusGLM<-glm(SuillusUn$Abundance ~ SuillusUn$TimeSinceFire + (1/Plot), data = SuillusUn)
SuillusAOV<-aov(SuillusGLM);summary(SuillusAOV)#0.285
SuillusTK<-TukeyHSD(SuillusAOV); SuillusTK #

ThelephoraGLM<-glm(ThelephoraUn$Abundance ~ ThelephoraUn$TimeSinceFire + (1/Plot), data = ThelephoraUn)
ThelephoraAOV<-aov(ThelephoraGLM);summary(ThelephoraAOV)#0.0163 *
ThelephoraTK<-TukeyHSD(ThelephoraAOV); ThelephoraTK #only 5-2

TomentellaGLM<-glm(TomentellaUn$Abundance ~ TomentellaUn$TimeSinceFire + (1/Plot), data = TomentellaUn)
TomentellaAOV<-aov(TomentellaGLM);summary(TomentellaAOV) #0.0447 *
TomentellaTK<-TukeyHSD(TomentellaAOV); TomentellaTK#5-3

TricholomaGLM<-glm(TricholomaUn$Abundance ~ TricholomaUn$TimeSinceFire + (1/Plot), data = TricholomaUn)
TricholomaAOV<-aov(TricholomaGLM);summary(TricholomaAOV) #0.00238 **
TricholomaTK<-TukeyHSD(TricholomaAOV); TricholomaTK# All years in comparison to 11

TuberGLM<-glm(TuberUn$Abundance ~ TuberUn$TimeSinceFire + (1/Plot), data = TuberUn)
TuberAOV<-aov(TuberGLM);summary(TuberAOV) #0.000407 ***
TuberTK<-TukeyHSD(TuberAOV); TuberTK #All in comparison to 5

TylosporaGLM<-glm(TylosporaUn$Abundance ~ TylosporaUn$TimeSinceFire + (1/Plot), data = TylosporaUn)
TylosporaAOV<-aov(TylosporaGLM);summary(TylosporaAOV) #0.00624 **
TylosporaTK<-TukeyHSD(TylosporaAOV); TylosporaTK #All in comparison to 5-2;11-2

WilcoxinaGLM<-glm(WilcoxinaUn$Abundance ~ WilcoxinaUn$TimeSinceFire + (1/Plot), data = WilcoxinaUn)
WilcoxinaAOV<-aov(WilcoxinaGLM);summary(WilcoxinaAOV) #9.88e-10 ***
WilcoxinaTK<-TukeyHSD(WilcoxinaAOV); WilcoxinaTK #All in comparison to 11



#Create a table and export the results-------------------------------------------------------------------------------------------------------
dir.create("Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn")

capture.output(summary(AmphinemaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Cennococum-glm.csv")
capture.output(summary(CortinariusAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Cortinarius-glm.csv")
capture.output(summary(EntolomaAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Entoloma-glm.csv")
capture.output(summary(HelvellosebacinaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/TSFunburn/Tables/Helvellosebacina-glm.csv")
capture.output(summary(HygrophorusAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Hygrophorus-glm.csv")
capture.output(summary(InocybeAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Inocybe-glm.csv")
capture.output(summary(PilodermaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Piloderma-glm.csv")
capture.output(summary(RhizopogonAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Rhizopogon-glm.csv")
capture.output(summary(RussulaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Russula-glm.csv")
capture.output(summary(SistotremaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Sistotrema-glm.csv")
capture.output(summary(SuillusAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Suillus-glm.csv")
capture.output(summary(ThelephoraAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Thelephora-glm.csv")
capture.output(summary(TomentellaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tomentella-glm.csv")
capture.output(summary(TricholomaAOV),file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tricholoma-glm.csv")
capture.output(summary(TuberAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tuber-glm.csv")
capture.output(summary(TylosporaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tylospora-glm.csv")
capture.output(summary(WilcoxinaAOV), file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Wilcoxina-glm.csv")

#Export Tukeys test------------------------------------------------------------------------------------------------------------------------
capture.output(AmphinemaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Amphinema-Tukey.csv")
capture.output(CortinariusTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Cortinarius-Tukey.csv")
capture.output(EntolomaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Entoloma-Tukey.csv")
capture.output(HelvellosebacinaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Helvellosebacina-Tukey.csv")
capture.output(HygrophorusTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Hygrophorus-Tukey.csv")
capture.output(InocybeTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Inocybe-Tukey.csv")
capture.output(PilodermaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Piloderma-Tukey.csv")
capture.output(RhizopogonTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Rhizopogon-Tukey.csv")
capture.output(RussulaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Russula-Tukey.csv")
capture.output(SistotremaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Sistotrema-Tukey.csv")
capture.output(SuillusTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Suillus-Tukey.csv")
capture.output(ThelephoraTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Thelephora-Tukey.csv")
capture.output(TomentellaTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tomentella-Tukey.csv")
capture.output(TricholomaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tricholoma-Tukey.csv")
capture.output(TuberTK,file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tuber-Tukey.csv")
capture.output(TylosporaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Tylospora-Tukey.csv")
capture.output(WilcoxinaTK, file="Analysis/RelativeAbundance/Ectomycorrhizal/Tables/TSFunburn/Wilcoxina-Tukey.csv")



























]
