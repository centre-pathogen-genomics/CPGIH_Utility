suppressMessages(suppressWarnings(library('tidyverse')))
suppressMessages(suppressWarnings(library('decontam')))
suppressMessages(suppressWarnings(library('qiime2R')))
suppressMessages(suppressWarnings(library('phyloseq')))

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]
metadata_file <- args[3]
decontam_column <- args[4]

phylo <- qza_to_phyloseq(features = input_file, metadata = metadata_file)

contamdf.prev <- decontam::isContaminant(phylo, method = "prevalence", neg = decontam_column)

noncontamOTUs <- contamdf.prev %>%
    filter(contaminant == FALSE) %>%
    row.names()

contamdf.prev %>%
    filter(contaminant == TRUE) %>%
    row.names() %>%
    write.csv('contaminant_otus.csv')

phylo_decontam <- prune_taxa(noncontamOTUs, phylo)

phylo_decontam %>%
    otu_table() %>%
    as.matrix() %>%
    biomformat::make_biom() %>%
    biomformat::write_biom(output_file)
