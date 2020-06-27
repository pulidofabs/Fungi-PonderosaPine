#!/usr/bin/bash
#SBATCH --nodes 4 --ntasks 12 --mem 36gb --time 12:00:00
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source activate qiime2-2019.10

awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' sh_refs_qiime_ver8_dynamic_s_02.02.2019_dev.fasta > sh_refs_qiime_ver8_dynamic_s_02.02.2019_dev_up.fasta

sed '67868 s/ //g' sh_refs_qiime_ver8_dynamic_s_02.02.2019_dev_up.fasta > sh_refs_qiime_ver8_dynamic_s_02.02.2019_dev_uppercase.fasta

qiime tools import \
--type FeatureData[Sequence] \
--input-path sh_refs_qiime_ver8_dynamic_s_02.02.2019_dev_uppercase.fasta \
--output-path unite-ver8-99-dynamic-seqs_s_02.02.2019.qza

qiime tools import \
--type FeatureData[Taxonomy] \
--input-path sh_taxonomy_qiime_ver8_dynamic_s_02.02.2019_dev.txt \
--output-path unite-ver8-99-tax-dynamic_s_02.02.2019.qza \
--input-format HeaderlessTSVTaxonomyFormat

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads unite-ver8-99-dynamic-seqs_s_02.02.2019.qza \
--i-reference-taxonomy unite-ver8-99-tax-dynamic_s_02.02.2019.qza \
--o-classifier classifier/unite-ver8-99-dynamic-classifier-02.02.2019.qza
