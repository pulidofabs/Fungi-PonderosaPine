#!/usr/bin/bash
#SBATCH --nodes 2 --ntasks 10 --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

qiime cutadapt trim-paired \
--i-demultiplexed-sequences Demux/Demux-paired.qza \
--p-adapter-f GCTGCGTTCTTCATCGATGC \
--p-adapter-f GCATCGATGAAGAACGCAGC \
--p-adapter-r CTTGGTCATTTAGAGGAAGTAA \
--p-adapter-r TTACTTCCTCTAAATGACCAAG \
--o-trimmed-sequences Demux/Demux-trimmed.qza \
--p-error-rate 0 \
--p-cores 10 \
--verbose

qiime demux summarize \
--i-data Demux/Demux-trimmed.qza \
--o-visualization Demux/Demux-trimmed.qzv
