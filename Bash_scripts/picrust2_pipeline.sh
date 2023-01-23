#!/bin/bash

############### ~~~~~~~~~~~ PICRUST2 Installation ~~~~~~~~~~~ ###############

# Before running this pipeline, you need to install and activate picrust2 conda environment
# To install --> conda install -c bioconda picrust2 
# Source --> https://anaconda.org/bioconda/picrust2
# To activate --> conda activate picrust2

############### ~~~~~~~~~~~ PICRUST2 PIPELINE ~~~~~~~~~~~ ###############

# Author: Zubair Khalid
# Date: 2021-06-01
# Description: This script will run the picrust2 pipeline
# The 

############### ~~~~~~~~~~~~~~Input files ~~~~~~~~~~~~~~ ################

# Input files required for the pipeline
# 1. Final Feature Table (converted in the qiime2-pipeline already using filtered-feature-table.qza)
# 2. Filtered Representative sequences 

# -----------------------------------------------------------------------

# Convert the representative sequences to a fasta file (Can also be downloaded from qiime2 view by using filtered-rep-seqs.qzv)
# The following commands need qiime2 to be activated. You can use another terminal to convert these. 
# Make sure the input directory name is correct and the files exist

ref_seqs_convert () {

    qiime tools export \
        --input-path inputs/filtered-rep-seqs.qza \
        --output-path inputs/rep-seqs
    
    if [ $? -eq 0 ]; then echo "Reference sequences converted to fasta file successfully"
        else "Error converting the reference sequences. Check if the directory name was correct and files exist"
    fi
}

ref_seqs_convert

# -----------------------------------------------------------------------

# The following is a full pipeline for PICRUST2
# The following commands need picrust2 to be activated. You can use another terminal to convert these.
# https://github.com/picrust/picrust2/wiki/Full-pipeline-script

picrust2_full_pipeline () {

    picrust2_pipeline.py \
        -s inputs/rep-seqs/dna-sequences.fasta \
        -i inputs/feature-table.biom \
        -o picrust2_outputs 

    if [ $? -eq 0 ]; then echo "PICRUST2 pipeline ran successfully"
        else "Error running the PICRUST2 pipeline"
    fi
}

picrust2_full_pipeline

# -----------------------------------------------------------------------

# The following commands will add descriptions to the output files

add_descriptions () {
    add_descriptions.py -i ./picrust2_outputs/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                        -o ./picrust2_outputs/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

    add_descriptions.py -i ./picrust2_outputs/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                        -o ./picrust2_outputs/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

    add_descriptions.py -i ./picrust2_outputs/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                        -o ./picrust2_outputs/pathways_out/path_abun_unstrat_descrip.tsv.gz
}

add_descriptions

# -----------------------------------------------------------------------
