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

noncontamOTUs <- rownames(contamdf.prev[contamdf.prev$contaminant == FALSE, ])

write.csv(rownames(contamdf.prev[contamdf.prev$contaminant == TRUE, ]), 
          file = "contaminant_otus.csv", 
          row.names = FALSE)

phylo_decontam <- prune_taxa(noncontamOTUs, phylo)

biom_obj <- biomformat::make_biom(as.matrix(otu_table(phylo_decontam)))
biomformat::write_biom(biom_obj, output_file)

