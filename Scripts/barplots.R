suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

# parse positional arguments
args <- commandArgs(trailingOnly = TRUE)

# make sure at least three arguments are provided
if (length(args) < 3) { stop('Error: The first three arguments are mandatory: namesOfSamplesToBePlotted, emuAbundanceProfile outputBarplot.pdf') }

# assign mandatory arguments
input_names <- args[1]
input_emu <- args[2]
output_plot <- args[3]

# assign optional arguments with default values if not provided
output_width <- if (length(args) >= 4) { as.numeric(args[4]) } else { 12 }
output_height <- if (length(args) >= 5) { as.numeric(args[5]) } else { 8 }

# import list of names to be plotted
names <- read_tsv(input_names, col_names = F)

# import EMu abundance file
df <- read_tsv(input_emu) %>%
  rename(Species = species)

# fill empty cells with zeroes
df[is.na(df)] <- 0

# pivot
df_long <- df %>%
  pivot_longer(cols = -Species, names_to = 'Sample', values_to = 'RA')

# filtering to only contain samples user wants plotted
df_long <- df_long %>%
    filter(Sample %in% names$X1)

# calculate mean RA of all species across samples to be plotted
meanRA <- df_long %>%
  group_by(Species) %>%
  summarise(MeanRA = mean(RA)) %>%
  arrange(-MeanRA)

# identify the 25 most abundant species by mean RA
topspecies <- meanRA %>%
  slice_max(order_by = MeanRA, n = 25)

# make colour palette for barplot
speciesColours <- c(hues::iwanthue(25), 'grey80')
names(speciesColours) <- c(sort(topspecies$Species), 'Other')

# plot
p <- df_long %>%
    mutate(Species = if_else(Species %in% topspecies$Species, Species, 'Other')) %>%
    mutate(Species = ordered(Species, levels = c(sort(topspecies$Species), 'Other'))) %>%
    mutate(Sample = ordered(Sample, levels = names$X1)) %>%
    group_by(Sample, Species) %>%
    summarise(RA= sum(RA)) %>%
    ggplot(aes(x = Sample, y = RA)) +
    geom_bar(aes(fill = Species), colour = 'black', stat = 'identity') +
    theme_bw() +
    scale_fill_manual(values = speciesColours) +
    labs(x = '',
         y = 'Relative Abundance') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

# save
pdf(output_plot, width = output_width, height = output_height)
p
dev.off()