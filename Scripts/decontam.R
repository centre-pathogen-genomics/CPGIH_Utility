# VERSION: one repo-wide number, stored in version.txt next to this script.
# Bump it (MAJOR = breaking change to inputs/outputs, MINOR = new feature, PATCH = bug fix) whenever any script changes
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] %in% c("--version", "-v")) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    version_file <- file.path(dirname(normalizePath(script_file)), "version.txt")
    version <- if (file.exists(version_file)) trimws(readLines(version_file, n = 1)) else "unknown (version.txt not found)"
    cat(basename(script_file), version, "\n")
    quit(status = 0)
}

suppressMessages(suppressWarnings(library('decontam')))
suppressMessages(suppressWarnings(library('qiime2R')))
suppressMessages(suppressWarnings(library('phyloseq')))

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

