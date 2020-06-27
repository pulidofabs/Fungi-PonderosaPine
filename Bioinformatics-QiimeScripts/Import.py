#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 4 --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

mkdir Demux #make a folder for where to store information

#Import files in format used by qiime
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path casava-18-paired-end-demultiplexed \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path Demux/Demux-paired.qza

qiime demux summarize \
--i-data Demux/Demux-paired.qza \
--o-visualization Demux/Demux-paired.qzv  

