#!/bin/bash

#  ================================== QIIME2 Pipeline =====================================  #

#  This script is to run the qiime2 pipeline for the microbiome analysis of 16S amplicons #
#  The script should be run in the directory containing the fastq.gz files and metadata.txt #
#  If your sequences are in folders, run extract_files.sh
#  The script is designed to prompt for manual inputs of several parameters such as ...
#       directory containing sequences, trim lengths, truncation lengths, seq. depths etc.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#  ====================================================================================  #
# ---- Author: Zubair Khalid  ---------------------------------------------------------- # 
# ---- Date: 30/12/2022  --------------------------------------------------------------- #
# ---- Version: 1.0  ------------------------------------------------------------------- # 
# ---- Source: https://docs.qiime2.org/2020.8/tutorials -------------------------------- #  
#  ====================================================================================  #     

# Load the module from Alabama Supercomputer if you're running the script there

# EDIT your USER here and add it to .bashrc.local file
# This step is not necessary

#source /home/MY_USER_NAME/.bashrc
#source /apps/profiles/modules_dmc.sh.dyn

module load qiime2/2023.2 
source activate qiime2-2023.2 

# If a more recent module has to be loaded, use the following command or tab after moudle load qiime
module spider qiime

# Make directory to manually dump the outputs from all the functions below eventually

mkdir -p qiime_outputs

# To get errors if a function fails
set -e

# Define the import data function
import_data() {
    
    echo "Enter the name of the directory containing the fastq.gz files"
    read fastq_dir
        
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path $fastq_dir \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt \
        --output-path demux-paired-end.qza
        
    if [ $? -eq 0 ]; then echo "Data import to qza completed successfully"
        else "Error importing the data. Check if the directory name was correct and files exist"
    fi
}

import_data  # calling the function

# Define the function for obtaining demultiplexing summary

demux_summarize() {
    
    qiime demux summarize \
        --i-data demux-paired-end.qza \
        --o-visualization demux_summary.qzv
        
    if [ $? -eq 0 ]; then
        echo "Demultiplexing summary is now ready to be viewed with QIIME2 View"
    fi
    
}

demux_summarize

# To denoise data, the demultiplexing summary will be analyzed and the inputs for denoising will be provided here
denoise() {
    echo "Please examine the denoising summary and input the forward and reverse trims."
        
    echo "Please input the length of forward primer (usually 20)"
    read trim_left_f
        
    echo "Please input the length of reverse primer (usually 20)"
    read trim_left_r
        
    echo "Enter the value for truncating length forward reads (240 for my file)"
    read trunc_len_f
        
    echo "Enter the value for truncating length for reverse reads (210 for my file)"
    read trunc_len_r
        
    qiime dada2 denoise-paired  \
        --i-demultiplexed-seqs demux-paired-end.qza \
        --p-trim-left-f $trim_left_f   --p-trim-left-r $trim_left_r \
        --p-trunc-len-f $trunc_len_f   --p-trunc-len-r $trunc_len_r  \
        --o-table table.qza   \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza
        
    if [ $? -eq 0 ]; then
        echo "Denoising completed successfully"
    fi
}

denoise

# Function to filter sequences and remove one ID that had only 5 feature counts compared to median of ~27000 for all others
filter_freq() {
    
    qiime feature-table filter-samples \
        --i-table table.qza \
        --p-min-frequency 1500 \
        --o-filtered-table sample-frequency-filtered-table.qza
        
    if [ $? -eq 0 ]; then
        echo "Samples with lower frequency filtered successfully"
    fi
}

filter_freq

# Function to summarize feature table
feature_table() {

    qiime feature-table summarize \
        --i-table sample-frequency-filtered-table.qza \
        --o-visualization table.qzv \
        --m-sample-metadata-file metadata.txt
        
    if [ $? -eq 0 ]; then
        echo "Feature table summarized successfully"
    fi
    
}

feature_table

#---> By visualizing feature table, I confirmed that there are 114 samples instead of 115 because of the low feature count

# Function to visualize sequences and save them as fasta from qiime view
rep_seqs_vis() {

    qiime feature-table tabulate-seqs \
        --i-data rep-seqs.qza \
        --o-visualization rep-seqs.qzv
    if [ $? -eq 0 ]; then
        echo "Representative sequences are ready to be visualized"
    fi   
}

rep_seqs_vis

# Tree for phylogenetic diversity

phylogeny() {

    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences rep-seqs.qza \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
    if [ $? -eq 0 ]; then
        echo "Rooted and unrooted trees has been generated"
    fi 
}

phylogeny


# Function to classify the feature tables based on taxa 
train_classifier() {

    qiime feature-classifier classify-sklearn \
        --i-classifier silva-138-99-nb-classifier.qza \
        --i-reads rep-seqs.qza  \
        --o-classification taxonomy.qza

    qiime metadata tabulate \
        --m-input-file taxonomy.qza   \
        --o-visualization taxonomy.qzv

    if [ $? -eq 0 ]; then
        echo "Taxonomy artifact and visualization was generated successfully"
    fi
}

train_classifier

# Filter the taxa to get only the bacterial sequences
filter_seqs () {

    qiime taxa filter-seqs \
        --i-sequences rep-seqs.qza \
        --i-taxonomy taxonomy.qza \
        --p-include d__,p__ \
        --p-exclude Archaea,Eukaryota,Unassigned,mitochondria,chloroplast \
        --o-filtered-sequences filtered-rep-seqs.qza
    
    if [ $? -eq 0 ]; then
        echo "Sequences containing Archaea, Eukaryota, Unassigned, mitochondria and chloroplasts have been filtered"
    fi
}

filter_seqs

# Function to visualize sequences and save them as fasta from qiime view
filtered_seqs_vis() {

    qiime feature-table tabulate-seqs \
        --i-data filtered-rep-seqs.qza \
        --o-visualization filtered-rep-seqs.qzv
    if [ $? -eq 0 ]; then
        echo "Filtered Representative sequences are ready to be visualized"
    fi   
}

filtered_seqs_vis

# Function to create taxa bar plots 
raw_taxa_barplot () {
    qiime taxa barplot  \
        --i-table sample-frequency-filtered-table.qza    \
        --i-taxonomy taxonomy.qza   \
        --m-metadata-file metadata.txt  \
        --o-visualization raw_taxa-bar-plots.qzv
    
    if [ $? -eq 0 ]; then
        echo "Barplot for taxonomic abundance is ready to be viewed"
    fi
}
raw_taxa_barplot

# Filter feature table based on taxa

final_filtered_table () {

    qiime taxa filter-table \
        --i-table sample-frequency-filtered-table.qza \
        --i-taxonomy taxonomy.qza \
        --p-include d__,p__ \
        --p-exclude Archaea,Eukaryota,Unassigned,mitochondria,chloroplast \
        --o-filtered-table final-feature-table.qza
    
    if [ $? -eq 0 ]; then
        echo "Features containing Archaea, Eukaryota, Unassigned, mitochondria and chloroplasts have been filtered"
    fi
}

final_filtered_table

id_renaming() {

    qiime feature-table rename-ids  \
        --i-table final-feature-table.qza   \
        --m-metadata-file metadata.txt  \
        --m-metadata-column label   \
        --p-strict True \
        --o-renamed-table final-feature-table.qza

    if [ $? -eq 0 ]; then
        echo "Feature table grouped by treatment"
    fi
}

id_renaming

# A new visualization of feature table with renamed ids and filtered features is warranted
final_table_vis() {

    qiime feature-table summarize \
        --i-table final-feature-table.qza \
        --o-visualization final-feature-table.qzv \
        --m-sample-metadata-file metadata.txt
        
    if [ $? -eq 0 ]; then
        echo "Feature table summarized successfully"
    fi
    
}

final_table_vis



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  Since the ids have been renamed, we should manually remove the column containing old IDs #
#############################################################################################

# Grouping by treatment and day by mean abundance

table_treat_grouping () {

    qiime feature-table group   \
        --i-table final-feature-table.qza   \
        --p-axis sample  \
        --m-metadata-file metadata.txt  \
        --m-metadata-column treatment   \
        --p-mode mean-ceiling   \
        --o-grouped-table treatment-grouped-table.qza
    
    if [ $? -eq 0 ]; then
        echo "Feature table grouped by treatment"
    fi

}

table_treat_grouping


table_dpi_grouping () {

    qiime feature-table group   \
        --i-table final-feature-table.qza   \
        --p-axis sample \
        --m-metadata-file metadata.txt  \
        --m-metadata-column sampling_day    \
        --p-mode mean-ceiling   \
        --o-grouped-table dpi-grouped-table.qza
    
    if [ $? -eq 0 ]; then
        echo "Feature table grouped by dpi"
    fi
}

table_dpi_grouping

# Function to create taxa bar plots 
filtered_taxa_barplot () {

    qiime taxa barplot  \
        --i-table final-feature-table.qza    \
        --i-taxonomy taxonomy.qza   \
        --m-metadata-file metadata.txt  \
        --o-visualization final-taxa-barplots.qzv
    
    if [ $? -eq 0 ]; then
        echo "Barplot for taxonomic abundance is ready to be viewed"
    fi
}

filtered_taxa_barplot

# Relative frequency of each feature per sample instead of absolute abundance

relative_freq () {

    qiime feature-table relative-frequency \
        --i-table final-feature-table.qza   \
        --o-relative-frequency-table rf-table.qza
    
    qiime feature-table summarize \
        --i-table rf-table.qza \
        --o-visualization rf-table.qzv \
        --m-sample-metadata-file metadata.txt

    if [ $? -eq 0 ]; then
        echo "Relative frequency table was created successfully"
    fi
}

relative_freq

# Alpha rarefaction plotting

alpha_rarefaction () {

    qiime diversity alpha-rarefaction   \
        --i-table final-feature-table.qza   \
        --i-phylogeny rooted-tree.qza   \
        --p-max-depth 20000 \
        --m-metadata-file metadata.txt  \
        --o-visualization alpha-rarefaction.qzv
    
    if [ $? -eq 0 ]; then
        echo "Relative frequency table was created successfully"
    fi
}

alpha_rarefaction

# Longitudinal analysis 
longitudinal () {

    qiime longitudinal volatility   \
        --i-table rf-table.qza   \
        --m-metadata-file metadata.txt  \
        --p-state-column numeric_dpi    \
        --p-default-group-column treatment  \
        --p-yscale log  \
        --o-visualization volatility.qzv    \
        --output-dir longitudinal
    
    if [ $? -eq 0 ]; then
        echo "Volatility file was created successfully"
    fi
}

longitudinal

# Yes, we are making collapsed tables but I am too tired to make the following lines of code pretty

collapsed_feature_tables () {

    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 5   --o-collapsed-table family-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 6   --o-collapsed-table genus-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 7   --o-collapsed-table species-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 1   --o-collapsed-table kingdom-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 2   --o-collapsed-table phylum-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 3   --o-collapsed-table class-table.qza 
    qiime taxa collapse   --i-table final-feature-table.qza   --i-taxonomy taxonomy.qza   --p-level 4   --o-collapsed-table order-table.qza

}

collapsed_feature_tables

# Relative frequency tables collapsed by taxa
collapsed_rel_freq_tables () {

    qiime feature-table relative-frequency --i-table genus-table.qza --o-relative-frequency-table rel-genus-table.qza
    qiime feature-table relative-frequency --i-table species-table.qza --o-relative-frequency-table rel-species-table.qza
    qiime feature-table relative-frequency --i-table order-table.qza --o-relative-frequency-table rel-order-table.qza
    qiime feature-table relative-frequency --i-table class-table.qza --o-relative-frequency-table rel-class-table.qza
    qiime feature-table relative-frequency --i-table phylum-table.qza --o-relative-frequency-table rel-phylum-table.qza
    qiime feature-table relative-frequency --i-table kingdom-table.qza --o-relative-frequency-table rel-kingdom-table.qza
    qiime feature-table relative-frequency --i-table family-table.qza --o-relative-frequency-table rel-family-table.qza
}

collapsed_rel_freq_tables

# Diversity analyses begin here!

diversity_core() {
    
    echo "Please specify the sampling depth after examining table.qzv in qiime view"
    read sampling_depth

    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table final-feature-table.qza \
        --p-sampling-depth $sampling_depth \
        --m-metadata-file metadata.txt \
        --output-dir diversity_core # A new directory will contain the outputs
        
    if [ $? -eq 0 ]; then 
        echo "Check the diversity core folder to see outputs"
    fi

}

diversity_core

alpha_group_sig () {

    qiime diversity alpha-group-significance \
        --i-alpha-diversity diversity_core/faith_pd_vector.qza \
        --m-metadata-file metadata.txt \
        --o-visualization diversity_core/faith-pd-group-significance.qzv

    qiime diversity alpha-group-significance \
        --i-alpha-diversity diversity_core/evenness_vector.qza \
        --m-metadata-file metadata.txt \
        --o-visualization diversity_core/evenness-group-significance.qzv


    if [ $? -eq 0 ]; then 
        echo "Faith's PD and Pileou's evenness artifacts generated successfully"
    fi
}

alpha_group_sig

beta_group_sig() {

    qiime diversity beta-group-significance \
        --i-distance-matrix diversity_core/unweighted_unifrac_distance_matrix.qza \
        --m-metadata-file metadata.txt \
        --m-metadata-column sampling_day \
        --o-visualization diversity_core/unweighted-unifrac-sampling-day-significance.qzv \
        --p-pairwise

    qiime diversity beta-group-significance \
        --i-distance-matrix diversity_core/unweighted_unifrac_distance_matrix.qza \
        --m-metadata-file metadata.txt \
        --m-metadata-column treatment \
        --o-visualization diversity_core/unweighted-unifrac-treatment-group-significance.qzv \
        --p-pairwise
    
        qiime diversity beta-group-significance \
        --i-distance-matrix diversity_core/weighted_unifrac_distance_matrix.qza \
        --m-metadata-file metadata.txt \
        --m-metadata-column sampling_day \
        --o-visualization diversity_core/weighted-unifrac-sampling-day-significance.qzv \
        --p-pairwise

    qiime diversity beta-group-significance \
        --i-distance-matrix diversity_core/weighted_unifrac_distance_matrix.qza \
        --m-metadata-file metadata.txt \
        --m-metadata-column treatment \
        --o-visualization diversity_core/weighted-unifrac-treatment-group-significance.qzv \
        --p-pairwise

    if [ $? -eq 0 ]; then 
        echo "Weighted and unweighted group significance artifacts generated successfully"
    fi
}

beta_group_sig

emperor_plots () {

    qiime emperor plot \
        --i-pcoa diversity_core/unweighted_unifrac_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --o-visualization diversity_core/unweighted-unifrac-emperor.qzv
    
    qiime emperor plot \
        --i-pcoa diversity_core/unweighted_unifrac_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --p-custom-axes numeric_dpi \
        --o-visualization diversity_core/unweighted-unifrac-dpi.qzv

    qiime emperor plot \
        --i-pcoa diversity_core/weighted_unifrac_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --o-visualization diversity_core/weighted-unifrac-emperor.qzv
    
    qiime emperor plot \
        --i-pcoa diversity_core/weighted_unifrac_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --p-custom-axes numeric_dpi \
        --o-visualization diversity_core/weighted-unifrac-dpi.qzv


    qiime emperor plot \
        --i-pcoa diversity_core/bray_curtis_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --o-visualization diversity_core/bray-curtis_pcoa.qzv

    qiime emperor plot \
        --i-pcoa diversity_core/bray_curtis_pcoa_results.qza \
        --m-metadata-file metadata.txt \
        --p-custom-axes numeric_dpi \
        --o-visualization diversity_core/bray-curtis_pcoa_dpi.qzv

    if [ $? -eq 0 ]; then 
        echo "Emperor plots generated successfully"
    fi
}

emperor_plots

volatility () {

    mkdir longitudinal

    qiime longitudinal volatility   \
        --i-table rel-genus-table.qza   \
        --m-metadata-file metadata.txt  \
        --m-metadata-file ./diversity_core/shannon_vector.qza   \
        --p-state-column numeric_dpi    \
        --p-default-group-column treatment  \
        --p-default-metric shannon_entropy  \
        --o-visualization ./longitudinal/volatility_shannon.qzv    \
    
    if [ $? -eq 0 ]; then
        echo "Relative frequency table was created successfully"
    fi
}

volatility

# Extract feature table in BIOM format. Not sure if I should use the collapsed table 

extract_biom () {

    mkdir -p extracted-feature-table

    qiime tools extract \
        --input-path final-feature-table.qza \
        --output-path extracted-feature-table
    if [ $? -eq 0 ]; then
        echo "Feature table was exported in biom format successfully"
    fi
}

extract_biom
