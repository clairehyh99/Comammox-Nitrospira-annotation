This pipeline provides a complete, end-to-end QIIME 2 workflow for the analysis of comammox Nitrospira amoA gene amplicons on the Spartan HPC environment.
It automates all steps from raw paired-end FASTQ files to quality-filtered, denoised ASVs and custom database–based taxonomy assignments, enabling reproducible, high-resolution annotation of comammox Nitrospira sequences.
The goal of this workflow is to produce amplicon sequence variants (ASVs) for comammox Nitrospira amoA genes with high accuracy and taxonomic precision, following an optimized version of the workflow published in:
Feng M., Lin Y., et al. (2025). “An optimized workflow for accurate taxonomic annotation of high-throughput comammox Nitrospira sequences.” (update the doi etc later)

Execution Environment

This pipeline is implemented as a single Slurm batch script compatible with the Spartan HPC system.

All stages (pre-processing and QIIME 2 processing) are run sequentially within one job.

Software modules

| Tool     | Version | Purpose                                      |
|-----------|----------|----------------------------------------------|
| **fastp** | 0.23.2  | Read quality control and trimming             |
| **FLASH** | 2.2.00  | Paired-end merging                            |
| **QIIME2** | 2022.11 | Denoising, taxonomic annotation, visualization |
| **GCC**   | 11.3.0  | Base compiler module required by Spartan HPC  |

Workflow Summary

1. Quality Filtering and Read Merging

a. Raw paired-end reads are processed using fastp v0.23.2 (Chen et al., 2018) and merged with FLASH v2.2.00 (Magoč et al., 2011) under the following criteria:

b. Reads are truncated where the average quality score < 20 within a 50 bp sliding window.

c. Truncated reads shorter than 50 bp and reads containing ambiguous bases (N) are discarded.

d. Overlapping sequences longer than 10 bp are merged, allowing a maximum mismatch ratio of 0.2.

e. Unmerged reads are discarded.

Each sample produces a *.merged.fastq.gz file. All operations are parallelized via Slurm (--cpus-per-task 16) with logs written to /path/to/processed.

2. Data Import and Quality Assessment

a. Merged single-end FASTQ files are automatically detected and converted into a QIIME 2 manifest.

b. Sequences are imported as SampleData[SequencesWithQuality] and summarized (demux summarize) to generate per-sample quality profiles.

3. Denoising and ASV Generation (DADA2)

Merged reads are denoised using the DADA2 algorithm within QIIME 2 to correct sequencing errors and infer unique ASVs. outputs include Representative ASV sequences (rep-seqs.qza), Feature table of ASV abundances (table.qza) and Denoising statistics (stats-dada2.qza)

4. Filtering of Low-Abundance ASVs

ASVs with fewer than 10 total reads are removed to eliminate sequencing noise and rare artifacts.

Filtered representative sequences (filtered-rep-seqs.qza) and tables (filtered_table.qza) are retained for downstream analysis.

5. Custom Database Annotation

User-defined reference FASTA and taxonomy files are imported into QIIME2:

custom_db.fasta → FeatureData[Sequence]

custom_taxonomy.txt → FeatureData[Taxonomy] (TSV Taxonomy Format)

6. Taxonomic Classification

Two complementary methods are used:

VSEARCH consensus classifier
Identity threshold: 0.8
Consensus threshold: 0.7
Outputs: taxonomy_vsearch.qza, vsearch_hits.qza

BLAST consensus classifier
Same parameters for independent validation
Outputs: taxonomy_blast.qza, blast_hits.qza

Both methods produce exported taxonomic tables and hit summaries for downstream curation.

7. Export outputs

The script automatically exports:

ASV abundance tables (feature-table.biom, TSV converted)

Representative sequences (FASTA)

DADA2 statistics and visualizations (.qzv)

VSEARCH / BLAST taxonomy results

Intermediate and final filtered feature tables

All outputs are stored under:
processed/01_fastp/
processed/02_merged/
processed/03_qiime2/
processed/logs/

Output summary

| Category          | Example Output                         | Description                                       |
|-------------------|-----------------------------------------|---------------------------------------------------|
| **Quality control** | `*.fastp.html`, `*.flash.log`          | Per-sample quality filtering and merging reports   |
| **QIIME2 data**     | `mjsample.qza`, `mjsample.qzv`         | Imported merged reads and quality summaries        |
| **DADA2 results**   | `rep-seqs.qza`, `table.qza`, `stats-dada2.qza` | Denoised ASVs, abundance table, and statistics  |
| **Filtered data**   | `filtered_table.qza`, `filtered-rep-seqs.qza` | ASVs retained after minimum read threshold filtering |
| **Taxonomy**        | `taxonomy_vsearch.qza`, `taxonomy_blast.qza` | Classifications via custom reference database     |
| **Exports**         | `asv_table.txt`, `filtered_asv_table.txt` | TSV-format abundance tables for downstream analysis |

Citation: 

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics. 2011 Nov 1;27(21):2957-63. doi: 10.1093/bioinformatics/btr507. Epub 2011 Sep 7. PMID: 21903629; PMCID: PMC3198573.

Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). https://doi.org/10.1038/s41587-019-0209-9
