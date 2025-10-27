This script was written by Mengmeng Feng as part of the publication: An optimized workflow for accurate taxonomic annotation of high-throughput comammox Nitrospira sequences.

The aim of this pipeline is to generate amplicon sequence variants (ASVs) of comammox (complete ammonia oxidizers) amoA genes, enabling a accurate，higher-resolution annotation of comammox Nitrospira . 

This workflow provides a complete end-to-end QIIME 2 pipeline. The raw data (paired-end sequences) were filtered using fastp v0.19.6 (Chen et al., 2018) and merged with FLASH v1.2.7 (Magoč et al., 2011) as follows. Reads were truncated at positions with an average quality score < 20 over a 50 bp sliding window; truncated reads < 50 bp and reads containing ambiguous characters were discarded. Overlapping sequences longer than 10 bp were assembled with a maximum mismatch ratio of 0.2.
