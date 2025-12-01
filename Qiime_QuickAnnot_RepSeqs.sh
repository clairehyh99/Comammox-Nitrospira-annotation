#!/usr/bin/env bash
set -e

###############################################
# Quick Annotation Script for Representative ASVs
# For users who only want to classify rep-seqs
# Author: Mengmeng Feng et al.
# Repository: Comammox-Nitrospira-annotation
###############################################

if [[ $# -lt 4 ]]; then
    echo "Usage: ./Qiime_QuickAnnot_RepSeqs.sh repseqs.qza custom_db.fasta taxonomy.txt output_dir"
    exit 1
fi

REPSEQ=$1
DB_FASTA=$2
DB_TAX=$3
OUT=$4

mkdir -p "$OUT"

echo "[1/3] Importing reference database..."
qiime tools import \
  --input-path "$DB_FASTA" \
  --output-path "$OUT/custom_db.qza" \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path "$DB_TAX" \
  --output-path "$OUT/custom_taxonomy.qza" \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat

echo "[2/3] Running VSEARCH classification..."
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$REPSEQ" \
  --i-reference-reads "$OUT/custom_db.qza" \
  --i-reference-taxonomy "$OUT/custom_taxonomy.qza" \
  --o-classification "$OUT/taxonomy_vsearch.qza" \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

echo "[3/3] Exporting results..."
qiime tools export \
  --input-path "$OUT/taxonomy_vsearch.qza" \
  --output-path "$OUT/taxonomy_exported"

echo "Quick representative-sequence annotation completed!"
echo "Results saved in: $OUT"
