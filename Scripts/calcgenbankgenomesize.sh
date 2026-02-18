
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

sed '1d' assembly_summary_genbank.txt | \
    cut -f 7,26 assembly_summary_genbank.txt | \
    csvtk summary -t -f genome_size:mean -f genome_size:median -f genome_size:stdev -g species_taxid | \
    taxonkit lineage -i 1 | \
    taxonkit reformat -i 5 | \
    csvtk cut -t -f 1,2,3,4,6 > taxid_genomelen_lineage.txt

rm -f assembly_summary_genbank.txt
