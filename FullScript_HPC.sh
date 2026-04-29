#!/usr/bin/env bash
#SBATCH --partition=YOUR_PARTITION   # replace with your HPC partition name
#SBATCH --nodes=1
#SBATCH --job-name=amplicon_qiime2
#SBATCH --account=YOUR_ACCOUNT       # replace with your HPC project account
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err

set -euo pipefail

# PATHS & PARAMETERS ## changed to your own directory
RAW_DIR="/path/to/rawData" #input dir
OUTROOT="/path/to/processed" #output dir
THREADS="${SLURM_CPUS_PER_TASK:-16}" #cpus

SAMPLE_TSV="/path/to/sample.tsv"                  
CUSTOM_DB_FASTA="/path/to/custom_db.fasta"         
CUSTOM_DB_TAX="/path/to/custom_taxonomy.txt"      

FP_DIR="${OUTROOT}/01_fastp"
MG_DIR="${OUTROOT}/02_merged"
Q2_DIR="${OUTROOT}/03_qiime2"
LOG_DIR="${OUTROOT}/logs"
mkdir -p "${FP_DIR}" "${MG_DIR}" "${Q2_DIR}" "${LOG_DIR}"

# load toolkits
module load GCC/11.3.0
module load fastp/0.23.2
module load FLASH/2.2.00
module load QIIME2/2022.11

echo "==> Listing samples from ${RAW_DIR}"
mapfile -t SAMPLES < <(find "${RAW_DIR}" -maxdepth 1 -type f -name "*.R1.raw.fastq.gz" -printf "%f\n" | awk -F '.' '{print $1}' | sort -V)
echo "==> Found ${#SAMPLES[@]} samples"

# fastp & cutadapt & FLASH
for s in "${SAMPLES[@]}"; do
  R1="${RAW_DIR}/${s}.CA209F_C576R.R1.raw.fastq.gz"
  R2="${RAW_DIR}/${s}.CA209F_C576R.R2.raw.fastq.gz"

  # fastp: Q20, 50bp sliding, min length 50, drop Ns
  fastp \
    -i "${R1}" -I "${R2}" \
    -o "${FP_DIR}/${s}.R1.fq.gz" -O "${FP_DIR}/${s}.R2.fq.gz" \
    --detect_adapter_for_pe \
    --cut_right --cut_window_size 50 --cut_mean_quality 20 \
    --length_required 50 \
    --n_base_limit 0 \
    --thread "${THREADS}" \
    --html "${LOG_DIR}/${s}.fastp.html" \
    --json "${LOG_DIR}/${s}.fastp.json" \
    > "${LOG_DIR}/${s}.fastp.log" 2>&1

  # Cutadapt: cutadapt fixed-length trimming
  cutadapt \
    -u 17 -U 19 \
    -o "cutadapt_output/${s}.R1.trimmed.fastq.gz" \
    -p "cutadapt_output/${s}.R2.trimmed.fastq.gz" \
    "${FP_DIR}/${s}.R1.fq.gz" \
    "${FP_DIR}/${s}.R2.fq.gz" \
    > "${LOG_DIR}/${s}.cutadapt.log" 2>&1

  # FLASH: min overlap 10, mismatch ratio 0.2
  flash \
    -m 10 \
    -x 0.20 \
    -t "${THREADS}" \
    -o "${s}" \
    -d "${MG_DIR}" \
    "cutadapt_output/${s}.R1.trimmed.fastq.gz" "cutadapt_output/${s}.R2.trimmed.fastq.gz" \
    > "${LOG_DIR}/${s}.flash.log" 2>&1

  if [[ -s "${MG_DIR}/${s}.extendedFrags.fastq" ]]; then
    gzip -f "${MG_DIR}/${s}.extendedFrags.fastq"
    mv "${MG_DIR}/${s}.extendedFrags.fastq.gz" "${MG_DIR}/${s}.merged.fastq.gz"
  else
    echo "[WARN] No merged reads for ${s}; see ${LOG_DIR}/${s}.flash.log" | tee -a "${LOG_DIR}/_warnings.txt"
  fi
done

echo "done. Merged FASTQs in: ${MG_DIR}"

# QIIME2: import + DADA2 + filters + custom DB taxonomy

# build a good list from merged files
MANIFEST="${Q2_DIR}/manifest_merged.csv"
echo "sample-id,absolute-filepath,direction" > "${MANIFEST}"
while IFS= read -r f; do
  bn="$(basename "${f}")"
  sid="${bn%%.*}"  
  echo "${sid},${f},forward" >> "${MANIFEST}"
done < <(find "${MG_DIR}" -maxdepth 1 -type f -name "*.merged.fastq.gz" | sort -V)

# Import merged reads
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --output-path "${Q2_DIR}/mjsample.qza" \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize quality
qiime demux summarize \
  --i-data "${Q2_DIR}/mjsample.qza" \
  --o-visualization "${Q2_DIR}/mjsample.qzv"

qiime tools export \
  --input-path "${Q2_DIR}/mjsample.qzv" \
  --output-path "${Q2_DIR}/mjsample_statistic"

# DADA2 (denoise-single on merged reads)
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${Q2_DIR}/mjsample.qza" \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences "${Q2_DIR}/rep-seqs.qza" \
  --o-table "${Q2_DIR}/table.qza" \
  --o-denoising-stats "${Q2_DIR}/stats-dada2.qza" \
  --p-n-threads "${THREADS}"

qiime metadata tabulate \
  --m-input-file "${Q2_DIR}/stats-dada2.qza" \
  --o-visualization "${Q2_DIR}/stats-dada2.qzv"

qiime tools export \
  --input-path "${Q2_DIR}/stats-dada2.qzv" \
  --output-path "${Q2_DIR}/stats2"

qiime feature-table summarize \
  --i-table "${Q2_DIR}/table.qza" \
  --o-visualization "${Q2_DIR}/table.qzv" \
  ${SAMPLE_TSV:+--m-sample-metadata-file "${SAMPLE_TSV}"}

qiime tools export \
  --input-path "${Q2_DIR}/table.qzv" \
  --output-path "${Q2_DIR}/table_stat"

# Export representative sequences
qiime tools export \
  --input-path "${Q2_DIR}/rep-seqs.qza" \
  --output-path "${Q2_DIR}/rep-seqs"

# Export ASV abundance (TSV)
mkdir -p "${Q2_DIR}/table"
qiime tools export \
  --input-path "${Q2_DIR}/table.qza" \
  --output-path "${Q2_DIR}/table"

biom convert \
  -i "${Q2_DIR}/table/feature-table.biom" \
  -o "${Q2_DIR}/asv_table.txt" \
  --table-type "OTU table" \
  --to-tsv

# Filter low-abundance ASVs (<10)
qiime feature-table filter-features \
  --i-table "${Q2_DIR}/table.qza" \
  --p-min-frequency 10 \
  --o-filtered-table "${Q2_DIR}/filtered_table.qza"

# Filter rep seqs to those retained
qiime feature-table filter-seqs \
  --i-data "${Q2_DIR}/rep-seqs.qza" \
  --i-table "${Q2_DIR}/filtered_table.qza" \
  --o-filtered-data "${Q2_DIR}/filtered-rep-seqs.qza"

# Visualize filtered table
qiime metadata tabulate \
  --m-input-file "${Q2_DIR}/filtered_table.qza" \
  --o-visualization "${Q2_DIR}/filtered_table.qzv"

qiime tools export \
  --input-path "${Q2_DIR}/filtered_table.qzv" \
  --output-path "${Q2_DIR}/filtered_table"

# Summarize filtered table (+metadata)
qiime feature-table summarize \
  --i-table "${Q2_DIR}/filtered_table.qza" \
  --o-visualization "${Q2_DIR}/filtered_table.qzv" \
  ${SAMPLE_TSV:+--m-sample-metadata-file "${SAMPLE_TSV}"}

qiime tools export \
  --input-path "${Q2_DIR}/filtered_table.qzv" \
  --output-path "${Q2_DIR}/filtered_table_stat"

# Export filtered rep seqs
qiime tools export \
  --input-path "${Q2_DIR}/filtered-rep-seqs.qza" \
  --output-path "${Q2_DIR}/filtered-rep-seqs"

# Export filtered abundance (TSV)
mkdir -p "${Q2_DIR}/filtered_table"
qiime tools export \
  --input-path "${Q2_DIR}/filtered_table.qza" \
  --output-path "${Q2_DIR}/filtered_table"

biom convert \
  -i "${Q2_DIR}/filtered_table/feature-table.biom" \
  -o "${Q2_DIR}/filtered_asv_table.txt" \
  --table-type "OTU table" \
  --to-tsv

# Custom DB taxonomy: vsearch & BLAST
# Import custom DB
qiime tools import \
  --input-path "${CUSTOM_DB_FASTA}" \
  --output-path "${Q2_DIR}/custom_db.qza" \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path "${CUSTOM_DB_TAX}" \
  --output-path "${Q2_DIR}/custom_taxonomy.qza" \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat

# VSEARCH consensus taxonomy
qiime feature-classifier classify-consensus-vsearch \
  --i-query "${Q2_DIR}/filtered-rep-seqs.qza" \
  --i-reference-reads "${Q2_DIR}/custom_db.qza" \
  --i-reference-taxonomy "${Q2_DIR}/custom_taxonomy.qza" \
  --o-classification "${Q2_DIR}/taxonomy_vsearch.qza" \
  --o-search-results "${Q2_DIR}/vsearch_hits.qza" \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

qiime tools export \
  --input-path "${Q2_DIR}/taxonomy_vsearch.qza" \
  --output-path "${Q2_DIR}/taxonomy_vsearch_exported"

qiime tools export \
  --input-path "${Q2_DIR}/vsearch_hits.qza" \
  --output-path "${Q2_DIR}/vsearch_hits_exported"

# BLAST consensus taxonomy
qiime feature-classifier classify-consensus-blast \
  --i-query "${Q2_DIR}/filtered-rep-seqs.qza" \
  --i-reference-reads "${Q2_DIR}/custom_db.qza" \
  --i-reference-taxonomy "${Q2_DIR}/custom_taxonomy.qza" \
  --o-classification "${Q2_DIR}/taxonomy_blast.qza" \
  --o-search-results "${Q2_DIR}/blast_hits.qza" \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

qiime tools export \
  --input-path "${Q2_DIR}/taxonomy_blast.qza" \
  --output-path "${Q2_DIR}/taxonomy_blast_exported"

echo "==> All done."
echo "    Preprocess: ${FP_DIR} / ${MG_DIR}"
echo "    QIIME2:     ${Q2_DIR}"
