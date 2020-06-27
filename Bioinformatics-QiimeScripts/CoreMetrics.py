#!/usr/bin/bash
#SBATCH --nodes 2 --ntasks 30 --mem 32gb --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

mkdir core-metrics-results

qiime diversity core-metrics-phylogenetic \
--i-phylogeny Phylo/Rooted-tree.qza \
--i-table DADA2/Table-dada2.qza \
--p-sampling-depth 10031 \
--m-metadata-file Metadata-PIPO.tsv \
--output-dir core-metrics-results

#Rarefy outside of core-metrics to see if it makes a difference----------------
qiime feature-table rarefy \
--i-table DADA2/Table-dada2.qza\
--p-sampling-depth 10031 \
--o-rarefied-table DADA2/Fungal-rarefied.qza

qiime feature-table summarize \
--i-table DADA2/Fungal-rarefied.qza \
--o-visualization DADA2/Fungal-rarefied.qzv \
--m-sample-metadata-file Metadata-PIPO.tsv



#CALCULATE SIMPSON DIVERSITY-------------------------------------------------
qiime diversity alpha \
--i-table core-metrics-results/rarefied_table.qza \
--p-metric simpson \
--o-alpha-diversity core-metrics-results/Simpson_vector.qza

qiime diversity alpha \
--i-table core-metrics-results/rarefied_table.qza \
--p-metric simpson_e \
--o-alpha-diversity core-metrics-results/SimpsonEvenness_vector.qza

#Simpson significance------------------------------------------------------------------------
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/Simpson_vector.qza \
--m-metadata-file  Metadata-PIPO.tsv \
--o-visualization core-metrics-results/Simpson_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/SimpsonEvenness_vector.qza \
--m-metadata-file  Metadata-PIPO.tsv \
--o-visualization core-metrics-results/SimpsonEvenness_vector.qzv

