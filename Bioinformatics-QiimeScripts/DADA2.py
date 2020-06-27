#!/usr/bin/bash
#SBATCH --nodes 2 --ntasks 30 --mem 32gb --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

mkdir DADA2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs Demux/Demux-trimmed.qza \
--p-trunc-len-f 150 \
--p-trunc-len-r 150 \
--p-trim-left-f 6 \
--p-trim-left-r 6 \
--p-trunc-q 0 \
--p-n-threads 30 \
--o-table DADA2/Table-dada2.qza \
--o-representative-sequences DADA2/RepSeqs-dada2.qza \
--o-denoising-stats DADA2/Stats-dada2.qza \
--verbose

qiime metadata tabulate \
--m-input-file DADA2/Stats-dada2.qza \
--o-visualization DADA2/Stats-dada2.qzv

qiime feature-table tabulate-seqs \
--i-data DADA2/RepSeqs-dada2.qza \
--o-visualization DADA2/RepSeqs-dada2.qzv

qiime feature-table summarize \
--i-table DADA2/Table-dada2.qza \
--o-visualization DADA2/Table-dada2.qzv \
