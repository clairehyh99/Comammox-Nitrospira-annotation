# QIIME 2 Workflow Script for Amplicon Data Processing and Custom Database Classification
# Author: Mengmeng Feng et al., 
# Date: [2025-10-18]

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
  --p-trunc-len 0 \
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
