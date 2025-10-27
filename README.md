This script was written by Mengmeng Feng as part of the publication: An optimized workflow for accurate taxonomic annotation of high-throughput comammox Nitrospira sequences.

The aim of this pipeline is to generate amplicon sequence variants (ASVs) of comammox (complete ammonia oxidizers) amoA genes, enabling a accurate，higher-resolution annotation of comammox Nitrospira . 

This workflow provides a complete end-to-end QIIME 2 pipeline. The raw data (paired-end sequences) were filtered using fastp v0.19.6 (Chen et al., 2018) and merged with FLASH v1.2.7 (Magoč et al., 2011) as follows. Reads were truncated at positions with an average quality score < 20 over a 50 bp sliding window; truncated reads < 50 bp and reads containing ambiguous characters were discarded. Overlapping sequences longer than 10 bp were assembled with a maximum mismatch ratio of 0.2.

The merged sequences were used in the pipeline for:
	1.	Data Import and Quality Assessment
Merged FASTQ files are imported into QIIME 2 and summarized to assess sequence quality. Quality reports are generated to guide the selection of truncation and filtering parameters for downstream analysis.
	2.	Denoising and ASV Generation (DADA2)
Raw reads are denoised using the DADA2 algorithm to correct sequencing errors and infer unique ASVs. Outputs include representative sequences, an ASV abundance table, and denoising statistics.
	3.	Filtering of Low-Abundance ASVs
ASVs with fewer than 10 reads are removed to eliminate sequencing noise and rare artifacts. The filtered representative sequences and abundance tables are then re-exported for subsequent analysis.
	4.	Annotation Using a Custom Database
A user-defined reference database (FASTA and taxonomy files) is imported and used to train a Naive Bayes classifier. The classifier is then applied to the filtered representative ASVs to assign taxonomic identities, generating an annotated taxonomy table.
	5.	Taxonomic Classification via VSEARCH
ASVs are additionally classified using the VSEARCH consensus-based method, employing user-defined identity (0.8) and consensus (0.7) thresholds. Both taxonomic assignments and sequence alignment results are exported for validation.
	6.	Taxonomic Classification via BLAST
A BLAST-based consensus classification is performed as a complementary validation step, using the same reference database and parameters as VSEARCH. Results are exported for downstream comparison and curation.

Please see more details in the publication!
