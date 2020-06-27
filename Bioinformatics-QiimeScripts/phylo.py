#!/usr/bin/bash
#SBATCH --nodes 2 --ntasks 30 --mem 32gb --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

mkdir Phylo

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences DADA2/RepSeqs-dada2.qza \
--o-alignment Phylo/Aligned-rep-seqs.qza \
--o-masked-alignment Phylo/Masked-aligned-rep-seqs.qza \
--o-tree Phylo/Unrooted-tree.qza \
--o-rooted-tree Phylo/Rooted-tree.qza

