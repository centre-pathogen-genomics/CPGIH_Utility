suppressMessages(suppressWarnings(library('tidyverse')))

args <- commandArgs(trailingOnly = TRUE)

input_directory <- args[1]

import_reports <- function(x){

    read_tsv(x, col_names = c('Species', 'Count')) %>%
        mutate(Sample = str_remove(x, '.*\\/')) %>%
        mutate(Sample = str_remove(Sample, '_report_species.tsv')) %>%
        mutate(Perc = 100*(Count/sum(Count))) %>%
        slice_max(Perc, n = 3, with_ties = F) %>%
        mutate(Name = paste0('Species', 1:3)) %>%
        mutate(Label = paste0(Species, ' (', Perc, '%)')) %>%
        pivot_wider(id_cols = Sample, names_from = Name, values_from = Label)

}

reports <- list.files(path = input_directory, pattern = '*_report_species.tsv', full.names = T)

map(reports, import_reports) %>% 
    bind_rows() %>%
    write_tsv(paste0(input_directory, 'top3species.tsv'))