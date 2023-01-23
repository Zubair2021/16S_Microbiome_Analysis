
qiime tools export --input-path qiime_outputs/genus-table.qza --output-path ./genus_table_biom
biom convert -i ./qiime_outputs/genus-table.biom -o feature-table.tsv --to-tsv
make_otu_network.py -i biom-feature-table/feature-table.biom -m metadata.txt -o otu_network
