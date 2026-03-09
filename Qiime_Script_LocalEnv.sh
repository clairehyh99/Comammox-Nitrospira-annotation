# QIIME 2 Workflow Script for Amplicon Data Processing and Custom Database Classification
# Author and affiliation anonymized for peer review
# Date: [2026-3-9]

# fastp: Quality-filtered 
mkdir -p fastp_output 
for sample in *_R1.fastq.gz; do
    sample_name=${sample%_R1.fastq.gz} 
    echo "Processing $sample_name with fastp..."

    fastp -i ${sample_name}_R1.fastq.gz -I ${sample_name}_R2.fastq.gz \
          -o fastp_output/${sample_name}_R1_clean.fastq.gz -O fastp_output/${sample_name}_R2_clean.fastq.gz \
          --cut_right --cut_right_window_size 50 --cut_right_mean_quality 20 \
          --length_required 50 --n_base_limit 0 \
          -h fastp_output/${sample_name}_fastp_report.html -j fastp_output/${sample_name}_fastp_report.json \
          -w 8
done

# cutadapt: Primer trimming (performed with trimming parameters corresponding to the actual primer lengths) 
mkdir -p cutadapt_output  
for sample in fastp_output/*_R1_clean.fastq.gz; do
    sample_name=$(basename "$sample" _R1_clean.fastq.gz)
    echo "Trimming fixed length from $sample_name with cutadapt..."

    cutadapt -u 17 -U 19 \
        -o cutadapt_output/${sample_name}_R1_trimmed.fastq.gz \
        -p cutadapt_output/${sample_name}_R2_trimmed.fastq.gz \
        fastp_output/${sample_name}_R1_clean.fastq.gz \
        fastp_output/${sample_name}_R2_clean.fastq.gz
done

# FLASH: Merged
mkdir -p flash_output  
for sample in cutadapt_output/*_R1_trimmed.fastq.gz; do
    sample_name=$(basename "$sample" _R1_trimmed.fastq.gz)
    echo "Merging $sample_name with FLASH..."

    flash cutadapt_output/${sample_name}_R1_trimmed.fastq.gz \
          cutadapt_output/${sample_name}_R2_trimmed.fastq.gz \
          -m 10 -x 0.2 \
          -d flash_output -o ${sample_name}

    gzip -c flash_output/${sample_name}.extendedFrags.fastq > flash_output/${sample_name}.extendedFrags.fastq.gz
    gzip -c flash_output/${sample_name}.notCombined_1.fastq > flash_output/${sample_name}.notCombined_1.fastq.gz
    gzip -c flash_output/${sample_name}.notCombined_2.fastq > flash_output/${sample_name}.notCombined_2.fastq.gz

done

# Import FASTQ data (Merged sequences)
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path mjsample.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize sequence quality
qiime demux summarize \
  --i-data mjsample.qza \
  --o-visualization mjsample.qzv

qiime tools export \
  --input-path mjsample.qzv \
  --output-path mjsample_statistic

# Denoising with DADA2

qiime dada2 denoise-single \
  --i-demultiplexed-seqs mjsample.qza \
  --p-trim-left 20 \
  --p-trunc-len 200 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza \
  --p-n-threads 8

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

qiime tools export \
  --input-path stats-dada2.qzv \
  --output-path stats2

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample.tsv

qiime tools export \
  --input-path table.qzv \
  --output-path table_stat

# Export representative sequences
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path rep-seqs

# Export ASV abundance table (as TSV)
qiime tools export \
  --input-path table.qza \
  --output-path table

biom convert \
  -i table/feature-table.biom \
  -o asv_table.txt \
  --table-type "OTU table" \
  --to-tsv

# Filter low-abundance ASVs (<10 reads)

qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 10 \
  --o-filtered-table filtered_table.qza

# Filter representative sequences based on the filtered feature table
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table filtered_table.qza \
  --o-filtered-data filtered-rep-seqs.qza

# Visualize filtered feature table (viewable via QIIME 2 View)
qiime metadata tabulate \
  --m-input-file filtered_table.qza \
  --o-visualization filtered_table.qzv

# Export .qzv visualization
qiime tools export \
  --input-path filtered_table.qzv \
  --output-path filtered_table

# Summarize ASV abundance and sample information
qiime feature-table summarize \
  --i-table filtered_table.qza \
  --o-visualization filtered_table.qzv \
  --m-sample-metadata-file sample.tsv

# Export summarized feature table
qiime tools export \
  --input-path filtered_table.qzv \
  --output-path filtered_table_stat

# Export filtered representative sequences
qiime tools export \
  --input-path filtered-rep-seqs.qza \
  --output-path filtered-rep-seqs

# Export filtered abundance table (as TSV)
qiime tools export \
  --input-path filtered_table.qza \
  --output-path filtered_table

biom convert \
  -i filtered_table/feature-table.biom \
  -o filtered_asv_table.txt \
  --table-type "OTU table" \
  --to-tsv



# Functional Gene Annotation Using a Custom Reference Database

qiime tools import \
  --input-path custom_db.fasta \
  --output-path custom_db.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path custom_taxonomy.txt \
  --output-path custom_taxonomy.qza \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat

# Taxonomic Classification Using Vsearch

qiime feature-classifier classify-consensus-vsearch \
  --i-query filtered-rep-seqs.qza \
  --i-reference-reads custom_db.qza \
  --i-reference-taxonomy custom_taxonomy.qza \
  --o-classification taxonomy_vsearch.qza \
  --o-search-results vsearch_hits.qza \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

qiime tools export \
  --input-path taxonomy_vsearch.qza \
  --output-path taxonomy_vsearch_exported

qiime tools export \
  --input-path vsearch_hits.qza \
  --output-path vsearch_hits_exported

# Taxonomic Classification Using BLAST

qiime feature-classifier classify-consensus-blast \
  --i-query filtered-rep-seqs.qza \
  --i-reference-reads custom_db.qza \
  --i-reference-taxonomy custom_taxonomy.qza \
  --o-classification taxonomy_blast.qza \
  --o-search-results blast_hits.qza \
  --p-perc-identity 0.8 \
  --p-min-consensus 0.7

qiime tools export \
  --input-path taxonomy_blast.qza \
  --output-path taxonomy_blast_exported
