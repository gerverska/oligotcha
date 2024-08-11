# Load libraries ####
library(ggplot2)
library(stringr)
library(dplyr)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Define input and output directories ####
in.path <- '02-align'
out <- '03-analyze'
plots <- file.path(out, 'plots')
rds <- file.path(out, 'rds')
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(plots, recursive = T)
dir.create(rds, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))

# Read in blastn results ####
columns <- c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
             'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

blockers <- read.delim(file.path(in.path, 'blockers.txt'),
                       col.names = columns) |> tax() |> filter(kingdom == 'Fungi')

# Extract the taxonomy associated with each entry in the database ####
ref <- file.path('data', '*silva-lsu-fungi.fa.gz')
paste("zgrep ';Fungi;'", ref, '>', file.path(out, 'headers.txt')) |> system()
ref.tax <- read.delim(file.path(out, 'headers.txt'),
                      col.names = 'stitle') |> tax()
ref.tax$source <- 'All fungi'

# Determine the number of fungal sequences in the reference database ####
total <- ref.tax |> nrow()

# Count the number of hits associated with each oligo ###
hits <- blockers |> group_by(qseqid) |> summarize(hits = n()) |>
    mutate(ra = hits / total)
hits$overlap <- hits$qseqid |> str_extract("[[:digit:]]+$") |> as.numeric()
best <- hits |> filter(hits == min(hits)) %>% .$qseqid

oligo <- blockers |> filter(qseqid == best)

# Subset the taxomnomy columns of the blastn output ####
oligo.tax <- oligo |> select(stitle, kingdom, phylum, class, order, family)
oligo.tax$source <- 'Blocking oligo BLAST hits'

# Combine both the reference and the output taxonomies ####
all.tax <- rbind(ref.tax, oligo.tax)

# Summarize counts at the phylum and class levels ####
phylum <- all.tax |> group_by(source, phylum) |> summarize(hits = n()) |>
    mutate(ra = hits / total) |> arrange(desc(ra)) |> 
    mutate(phylum = factor(phylum,
                          levels = unique(phylum)))

class <- all.tax |> group_by(source, class) |> summarize(hits = n()) |> 
    mutate(ra = hits / total) |> arrange(desc(ra)) |> 
    mutate(class = factor(class,
                          levels = unique(class)))

# Count the number of non-target hits of a given length and starting position
start.end <- oligo |> group_by(qstart, qend) |> summarize(hits = n()) |>
    mutate(ra = hits / total) |> filter(hits > 1)

# Plot the number of non-target hits for each phylum ####
overlap.plot <- ggplot(hits, aes(x = overlap, y = hits)) +
    geom_bar(stat = 'identity',
             fill = 'black',
             color = 'white') +
    xlab('\nBlocking oligo overlap with primer (bp)') +
    ylab('BLAST hits to database\n') +
    scale_x_continuous(n.breaks = 15) +
    scale_y_continuous(n.breaks = 15) +
    theme_classic() +
    base.theme
file.path(plots, 'overlap.png') |> ggsave(overlap.plot)
file.path(rds, 'overlap.rds') |> saveRDS(overlap.plot, file = _)

# Plot the number of non-target hits for each phylum ####
phylum.plot <- ggplot(phylum, aes(x = phylum, y = ra)) +
    geom_bar(stat = 'identity',
             fill = 'black',
             color = 'white') +
    facet_grid(rows = vars(source), scales = 'free_y') +
    xlab('\nPhylum') +
    ylab('Proportion of total database\n') +
    scale_y_continuous(n.breaks = 6,
                       trans = 'sqrt',
                       labels = function(x) format(x, scientific = FALSE)) +
    theme_classic() +
    base.theme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
file.path(plots, 'phylum.png') |> ggsave(phylum.plot)
file.path(rds, 'phylum.rds') |> saveRDS(phylum.plot, file = _)

# Plot the number of non-target hits for each class ####
class.plot <- ggplot(class, aes(x = class, y = ra)) +
    geom_bar(stat = 'identity',
             fill = 'black',
             color = 'white') +
    facet_grid(rows = vars(source), scales = 'free_y') +
    xlab('\nClass') +
    ylab('Proportion of total database\n') +
    scale_y_continuous(n.breaks = 12,
                       trans = 'sqrt',
                       labels = function(x) format(x, scientific = FALSE)) +
    theme_classic() +
    base.theme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
file.path(plots, 'class.png') |> ggsave(class.plot, height = 9.5)
file.path(rds, 'class.rds') |> saveRDS(class.plot, file = _)

# Plot the distribution of query alignment lengths ####
length.plot <- ggplot(oligo, aes(x = length)) +
    geom_histogram(binwidth = 1,
                   fill = 'black',
                   color = 'white') +
    xlab('\nAlignment length (bp)') +
    ylab('BLAST hits to database\n') +
    scale_x_continuous(n.breaks = 20) +
    scale_y_continuous(n.breaks = 15) +
    theme_classic() +
    base.theme
file.path(plots, 'length.png') |> ggsave(length.plot)
file.path(rds, 'length.rds') |> saveRDS(length.plot, file = _)

# Plot the distribution of alignment percent identity ####
pident.plot <- ggplot(oligo, aes(x = pident)) +
    geom_histogram(binwidth = 1,
                   fill = 'black',
                   color = 'white') +
    xlab("\nPercent identity") +
    ylab('BLAST hits to database\n') +
    scale_y_continuous(n.breaks = 15) +
    theme_classic() +
    base.theme
file.path(plots, 'pident.png') |> ggsave(pident.plot)
file.path(rds, 'pident.rds') |> saveRDS(pident.plot, file = _)

# Plot the relationship between query alignment start and end ####
start.end.plot <- ggplot(start.end,
                            aes(x = qstart, y = qend, color = hits)) +
    geom_point(size = 3) +
    scale_color_viridis_c() +
    xlab("\nQuery start position") +
    ylab("Query end position\n") +
    labs(color = 'BLAST hits') +
    scale_x_continuous(n.breaks = 20) +
    scale_y_continuous(n.breaks = 20) +
    theme_classic() +
    base.theme
file.path(plots, 'start-end.png') |> ggsave(start.end.plot, width = 8)
file.path(rds, 'start_end.rds') |> saveRDS(start.end.plot, file = _)
