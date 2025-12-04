# Comammox *Nitrospira* amoA Amplicon Analysis Pipeline  
### Local QIIME 2 Wrapper Script (FASTQ → DADA2 → Filter → Custom Annotation)

This repository provides a complete **local workstation pipeline** for high-resolution analysis of *comammox Nitrospira* **amoA** gene amplicons.

The pipeline automates all steps from **raw paired-end FASTQ** files to **denoised ASVs** and **custom database–based taxonomic annotation**, following an optimized workflow described in:

> **Feng M., Lin Y., et al. (2025)**.  
> *An optimized workflow for accurate taxonomic annotation of high-throughput comammox Nitrospira amoA sequences.*  
> *(DOI coming soon)*

---


# Overview

This pipeline is implemented as a single Bash wrapper:

Two modes are supported:

| Mode | Description |
|------|-------------|
| `full` | Run complete workflow from raw FASTQ → ASVs → annotation |
| `annot` | Only annotation on existing `filtered-rep-seqs.qza` |


---

# Execution Environment

Recommended conda setup:

```bash
conda create -n qiime2-amplicon-2024.5 python=3.10
conda activate qiime2-amplicon-2024.5
mamba install fastp cutadapt flash biom-format

Required software
| Tool            | Purpose                 | Link                                                                     |
| --------------- | ----------------------- | ------------------------------------------------------------------------ |
| **QIIME 2**     | ASV inference, taxonomy | [https://qiime2.org](https://qiime2.org)                                 |
| **fastp**       | FASTQ QC                | [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)   |
| **cutadapt**    | Primer trimming         | [https://cutadapt.readthedocs.io](https://cutadapt.readthedocs.io)       |
| **FLASH**       | Paired-end merging      | [https://ccb.jhu.edu/software/FLASH](https://ccb.jhu.edu/software/FLASH) |
| **biom-format** | OTU table conversion    | [https://biom-format.org](https://biom-format.org)                       |

Input requirements

project/
├── raw_fastq/
│   ├── F1.CA209F_C576R.R1.raw.fastq.gz
│   ├── F1.CA209F_C576R.R2.raw.fastq.gz
│   └── ...
├── comammox_db.fasta
├── comammox_tax.txt
├── sample_metadata.tsv   (optional)
└── Wrap_up_FASTQ-DADA2-FILTER-ANNOTATION.sh

FASTQ file naming pattern:
<sample>.CA209F_C576R.R1.raw.fastq.gz
<sample>.CA209F_C576R.R2.raw.fastq.gz

Custom DB format:

Feature ID    Taxon
CUS38424       A; A.1
MH444516       A; A.2

Usage
# Entire pipeline
bash Wrap_up_FASTQ-DADA2-FILTER-ANNOTATION.sh full \
  --raw-dir ./raw_fastq \
  --outdir ./out_amplicon \
  --sample-tsv sample_metadata.tsv \
  --custom-db-fasta comammox_db.fasta \
  --custom-db-tax comammox_tax.txt \
  --threads 8

# Annotation only
bash Wrap_up_FASTQ-DADA2-FILTER-ANNOTATION.sh annot \
  --repseq-qza out_amplicon/04_qiime2/filtered-rep-seqs.qza \
  --annot-outdir out_amplicon/04_qiime2/comammox_annotation \
  --custom-db-fasta comammox_db.fasta \
  --custom-db-tax comammox_tax.txt

Workflow Summary
1. fastp quality filtering (50 bp window, Q < 20)

Removes low-quality reads, N-containing reads, short fragments.

2. cutadapt trimming

Removes fixed-length primer sequence (17 bp / 19 bp).

3. FLASH merging

Minimum overlap: 10 bp

Mismatch ratio: 0.2

Output: *.merged.fastq.gz

4. QIIME2 import

Auto-generates manifest with absolute paths.

5. DADA2

Produces:

rep-seqs.qza

table.qza

dada2_stats.qza

6. Filter ASVs (min frequency ≥ 10)

Outputs:

filtered-rep-seqs.qza

filtered_table.qza

7. Custom comammox annotation

Two classifiers:

Method	Thresholds
VSEARCH	identity ≥ 0.8, consensus ≥ 0.7
BLAST	identity ≥ 0.8, consensus ≥ 0.7

Exports TSV tables for curating comammox clades.

📂 Output Structure
out_amplicon/
├── 01_fastp/
├── 02_cutadapt/
├── 03_merged/
└── 04_qiime2/
    ├── manifest.tsv
    ├── mjsample.qza
    ├── rep-seqs.qza
    ├── table.qza
    ├── dada2_stats.qza
    ├── filtered-rep-seqs.qza
    ├── filtered_table.qza
    ├── filtered_asv_table.txt
    └── comammox_annotation/
        ├── taxonomy_vsearch.qza
        ├── taxonomy_vsearch_exported/
        ├── taxonomy_blast.qza
        └── taxonomy_blast_exported/

Citation: 

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics. 2011 Nov 1;27(21):2957-63. doi: 10.1093/bioinformatics/btr507. Epub 2011 Sep 7. PMID: 21903629; PMCID: PMC3198573.

Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). https://doi.org/10.1038/s41587-019-0209-9
