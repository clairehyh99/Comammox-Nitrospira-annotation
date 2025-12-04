#!/usr/bin/env bash
set -euo pipefail
# Modes:
#   full  – FASTQ → DADA2 → filtered-rep-seqs → annotation
#   annot – annotation only
MODE="${1:-}"
if [[ -z "${MODE}" ]]; then
  echo "Usage:"
  echo "  $0 full  --raw-dir DIR --outdir DIR --sample-tsv FILE --custom-db-fasta FILE --custom-db-tax FILE [--threads N]"
  echo "  $0 annot --repseq-qza QZA --annot-outdir DIR --custom-db-fasta FILE --custom-db-tax FILE [--threads N]"
  exit 1
fi
shift
THREADS=8  
RAW_DIR=""
OUTROOT=""
SAMPLE_TSV=""
CUSTOM_DB_FASTA=""
CUSTOM_DB_TAX=""
REPSEQ_QZA=""
ANNOT_OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-dir) RAW_DIR="$2"; shift 2 ;;
    --outdir) OUTROOT="$2"; shift 2 ;;
    --sample-tsv) SAMPLE_TSV="$2"; shift 2 ;;
    --custom-db-fasta) CUSTOM_DB_FASTA="$2"; shift 2 ;;
    --custom-db-tax) CUSTOM_DB_TAX="$2"; shift 2 ;;
    --repseq-qza) REPSEQ_QZA="$2"; shift 2 ;;
    --annot-outdir) ANNOT_OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *)
      echo "[ERR] Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

########################## FUNCTION: Annotation module (cmx) ###########################################
run_annotation_only() {
  local repseq="$1"
  local db_fasta="$2"
  local db_tax="$3"
  local outdir="$4"
  mkdir -p "$outdir"
  echo "====================================================="
  echo "[ANNOT] Running Comammox Annotation"
  echo "====================================================="
  # Import DB
  qiime tools import \
    --input-path "${db_fasta}" \
    --output-path "${outdir}/custom_db.qza" \
    --type 'FeatureData[Sequence]'
  qiime tools import \
    --input-path "${db_tax}" \
    --output-path "${outdir}/custom_tax.qza" \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat
  # VSEARCH
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "${repseq}" \
    --i-reference-reads "${outdir}/custom_db.qza" \
    --i-reference-taxonomy "${outdir}/custom_tax.qza" \
    --p-perc-identity 0.8 \
    --p-min-consensus 0.7 \
    --o-classification "${outdir}/taxonomy_vsearch.qza" \
    --o-search-results "${outdir}/vsearch_hits.qza"
  qiime tools export \
    --input-path "${outdir}/taxonomy_vsearch.qza" \
    --output-path "${outdir}/taxonomy_vsearch_exported"
  # BLAST
  qiime feature-classifier classify-consensus-blast \
    --i-query "${repseq}" \
    --i-reference-reads "${outdir}/custom_db.qza" \
    --i-reference-taxonomy "${outdir}/custom_tax.qza" \
    --p-perc-identity 0.8 \
    --p-min-consensus 0.7 \
    --o-classification "${outdir}/taxonomy_blast.qza" \
    --o-search-results "${outdir}/blast_hits.qza"
  qiime tools export \
    --input-path "${outdir}/taxonomy_blast.qza" \
    --output-path "${outdir}/taxonomy_blast_exported"
  echo "====================================================="
  echo "[ANNOT] Completed. Output in: $outdir"
  echo "====================================================="
}


################################# FUNCTION: Full pipeline (FASTQ → ASV → annotation) #########################################
run_full_pipeline() {
  local raw_dir="$1"
  local outroot="$2"
  local sample_tsv="$3"
  local db_fasta="$4"
  local db_tax="$5"
  local threads="$6"
  local FP_DIR="${outroot}/01_fastp"
  local CUT_DIR="${outroot}/02_cutadapt"
  local MG_DIR="${outroot}/03_merged"
  local Q2_DIR="${outroot}/04_qiime2"
  mkdir -p "$FP_DIR" "$CUT_DIR" "$MG_DIR" "$Q2_DIR"
  echo "====================================================="
  echo "[FULL] Starting Full Pipeline"
  echo "====================================================="
  mapfile -t SAMPLES < <(
    find "$raw_dir" -name "*.R1.raw.fastq.gz" | sed 's/.CA209F_C576R.R1.raw.fastq.gz//'
  )
  for s in "${SAMPLES[@]}"; do
    base=$(basename "$s")
    R1="${s}.CA209F_C576R.R1.raw.fastq.gz"
    R2="${s}.CA209F_C576R.R2.raw.fastq.gz"
    echo "[FASTP] $base"
    fastp \
      -i "$R1" -I "$R2" \
      -o "$FP_DIR/${base}.R1.fq.gz" \
      -O "$FP_DIR/${base}.R2.fq.gz" \
      --detect_adapter_for_pe \
      --cut_right --cut_window_size 50 --cut_mean_quality 20 \
      --length_required 50 \
      --thread "$threads"
    echo "[CUTADAPT] $base"
    cutadapt \
      -u 17 -U 19 \
      -o "$CUT_DIR/${base}.R1.trimmed.fastq.gz" \
      -p "$CUT_DIR/${base}.R2.trimmed.fastq.gz" \
      "$FP_DIR/${base}.R1.fq.gz" \
      "$FP_DIR/${base}.R2.fq.gz"
    echo "[FLASH] $base"
    flash \
      "$CUT_DIR/${base}.R1.trimmed.fastq.gz" \
      "$CUT_DIR/${base}.R2.trimmed.fastq.gz" \
      -o "$base" \
      -d "$MG_DIR" \
      -m 10 -x 0.2 -t "$threads"
    if [[ -f "$MG_DIR/${base}.extendedFrags.fastq" ]]; then
      gzip -f "$MG_DIR/${base}.extendedFrags.fastq"
      mv "$MG_DIR/${base}.extendedFrags.fastq.gz" "$MG_DIR/${base}.merged.fastq.gz"
    fi
  done
  # Import to QIIME2
  echo "[QIIME2] Import and DADA2"
  MANIFEST="$Q2_DIR/manifest.tsv"
echo -e "sample-id\tabsolute-filepath" > "$MANIFEST"
 
for f in "$MG_DIR"/*.merged.fastq.gz; do
    bn=$(basename "$f")
    sid="${bn%%.*}"
    echo -e "$sid\t$(readlink -f "$f")" >> "$MANIFEST"
done
 
dos2unix "$MANIFEST"
 
  qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
  --input-path "$Q2_DIR/manifest.tsv" \
  --output-path "$Q2_DIR/mjsample.qza" \
  --input-format SingleEndFastqManifestPhred33V2
  qiime dada2 denoise-single \
    --i-demultiplexed-seqs "$Q2_DIR/mjsample.qza" \
    --p-trim-left 20 \
    --p-trunc-len 200 \
    --p-n-threads "$threads" \
    --o-representative-sequences "$Q2_DIR/rep-seqs.qza" \
    --o-table "$Q2_DIR/table.qza" \
    --o-denoising-stats "$Q2_DIR/dada2_stats.qza"
  # Filter low abundance
  qiime feature-table filter-features \
    --i-table "$Q2_DIR/table.qza" \
    --p-min-frequency 10 \
    --o-filtered-table "$Q2_DIR/filtered_table.qza"
  qiime feature-table filter-seqs \
    --i-data "$Q2_DIR/rep-seqs.qza" \
    --i-table "$Q2_DIR/filtered_table.qza" \
    --o-filtered-data "$Q2_DIR/filtered-rep-seqs.qza"

  ############################################################ Wrap-up: run annotation module###########################################################
  run_annotation_only \
    "$Q2_DIR/filtered-rep-seqs.qza" \
    "$db_fasta" \
    "$db_tax" \
    "$Q2_DIR/comammox_annotation"
  echo "====================================================="
  echo "[FULL] Completed."
  echo "Results in: $outroot"
  echo "====================================================="
}

############################################## Main logic
case "$MODE" in
  full)
    run_full_pipeline "$RAW_DIR" "$OUTROOT" "$SAMPLE_TSV" \
                      "$CUSTOM_DB_FASTA" "$CUSTOM_DB_TAX" "$THREADS"
    ;;
  annot)
    run_annotation_only "$REPSEQ_QZA" "$CUSTOM_DB_FASTA" "$CUSTOM_DB_TAX" "$ANNOT_OUTDIR"
    ;;
  *)
    echo "[ERR] Unknown MODE: $MODE"
    exit 1
    ;;
esac