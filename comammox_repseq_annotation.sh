#!/usr/bin/env bash
set -e

# =====================================================================
# One-Click Wrapper Script for Comammox Representative Sequence Annotation 
# Author names and affiliations have been removed for double-anonymous peer review. 
# Date: [2026-3-9]
#
# This script performs:
#   1. Importing the custom comammox amoA reference database
#   2. VSEARCH-based taxonomic classification
#   3. BLAST-based taxonomic classification
#   4. Exporting all annotation results
#
# Usage:
#   ./comammox_repseq_annotation.sh filtered-rep-seqs.qza custom_db.fasta custom_taxonomy.txt output_dir
# =====================================================================

# -------------------- Input Arguments --------------------
REPSEQ=$1          # filtered representative sequences (QZA)
DB_FASTA=$2        # custom_db.fasta
DB_TAX=$3          # custom_taxonomy.txt
OUTDIR=$4          # output directory

if [[ $# -lt 4 ]]; then
  echo "Usage: ./comammox_repseq_annotation.sh filtered-rep-seqs.qza custom_db.fasta custom_taxonomy.txt output_dir"
  exit 1
fi

mkdir -p "$OUTDIR"

echo "================================================================="
echo "  Running One-Click Comammox Annotation Pipeline"
echo "  Representative sequences : $REPSEQ"
echo "  Database FASTA          : $DB_FASTA"
echo "  Database Taxonomy       : $DB_TAX"
echo "  Output directory        : $OUTDIR"
echo "================================================================="


# -------------------- Step 1: Import Reference Database --------------------
echo "[1/3] Importing custom comammox amoA reference database..."

qiime tools import \
  --input-path "$DB_FASTA" \
  --output-path "$OUTDIR/custom_db.qza" \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path "$DB_TAX" \
  --output-path "$OUTDIR/custom_taxonomy.qza" \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat


# -------------------- Step 2: VSEARCH Classification --------------------
echo "[2/3] Running VSEARCH-based taxonomic classification..."

qiime feature-classifier classify-consensus-vsearch \
  --i-query "$REPSEQ" \
  --i-reference-reads "$OUTDIR/custom_db.qza" \
  --i-reference-taxonomy "$OUTDIR/custom_taxonomy.qza" \
  --o-classification "$OUTDIR/taxonomy_vsearch.qza" \
  --o-search-results "$OUTDIR/vsearch_hits.qza" \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

# Export VSEARCH results
qiime tools export \
  --input-path "$OUTDIR/taxonomy_vsearch.qza" \
  --output-path "$OUTDIR/taxonomy_vsearch_exported"

qiime tools export \
  --input-path "$OUTDIR/vsearch_hits.qza" \
  --output-path "$OUTDIR/vsearch_hits_exported"


# -------------------- Step 3: BLAST Classification --------------------
echo "[3/3] Running BLAST-based taxonomic classification..."

qiime feature-classifier classify-consensus-blast \
  --i-query "$REPSEQ" \
  --i-reference-reads "$OUTDIR/custom_db.qza" \
  --i-reference-taxonomy "$OUTDIR/custom_taxonomy.qza" \
  --o-classification "$OUTDIR/taxonomy_blast.qza" \
  --o-search-results "$OUTDIR/blast_hits.qza" \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

# Export BLAST results
qiime tools export \
  --input-path "$OUTDIR/taxonomy_blast.qza" \
  --output-path "$OUTDIR/taxonomy_blast_exported"


echo "================================================================="
echo "  Comammox Annotation Completed Successfully!"
echo "  Results saved to: $OUTDIR"
echo "================================================================="

