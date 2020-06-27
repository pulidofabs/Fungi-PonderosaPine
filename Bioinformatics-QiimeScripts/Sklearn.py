#!/usr/bin/bash
#SBATCH --nodes 2 --ntasks 30 --mem 32gb --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

mkdir Taxonomy

qiime feature-classifier classify-sklearn \
--i-classifier Classifier2020-Fun/unite-ver8-99-classifier-04.02.2020.qza \
--i-reads DADA2/RepSeqs-dada2.qza \
--p-read-orientation same \
--o-classification Taxonomy/Fungal-taxonomy-paired.qza \
--verbose

qiime metadata tabulate \
--m-input-file Taxonomy/Fungal-taxonomy-paired.qza \
--o-visualization Taxonomy/Fungal-taxonomy-paired.qzv


#Create taxonomy barplots
qiime taxa barplot \
--i-table DADA2/Fungal-Table-dada2-NS.qza \
--i-taxonomy Taxonomy/Fungal-taxonomy-paired.qza \
--m-metadata-file  Metadata-PIPO.tsv \
--o-visualization Taxonomy/Fungal-Tax-Barplot.qzv






