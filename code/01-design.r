# Accept an argument for the number of threads and check it for validity ####
threads <- commandArgs(T) |> as.integer()

if(is.na(threads) == T){
    stop('Argument was converted to NA')
}
if(length(threads) < 1){
    stop('Please specify the number of threads to launch')
}
if(length(threads) > 1){
    stop('Too many arguments have been provided')
}
if(is.numeric(threads) == F){
    stop('Only numeric arguments are accepted')
}
if(threads < 1){
    stop('At least one thread is needed')
} else {
    cat(threads, 'threads requested', '\n')
}

# Load packages ####
library(Biostrings)
library(dplyr)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Define input and output directories ####
in.path <- 'data'
out <- '01-design'
db <- file.path(out, 'db')
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(db, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))

# Make the host database ####
db.args <- paste(
    '-in', list.files(in.path, pattern = '*host.fa', full.names = T),
    '-dbtype nucl',
    '-title "Host"',
    '-out', file.path(db, 'host'),
    '-logfile', file.path(logs, 'host.txt')
)
system2('makeblastdb', args = db.args)

# Align the fungal primer to the host database ####
blastn.args <- paste(
    '-query', list.files(in.path, pattern = '*primer.fa', full.names = T),
    '-db', file.path(db, 'host'),
    '-dust no',
    '-evalue 1000',
    '-perc_identity 50',
    '-word_size 10',
    '-max_target_seqs 10',
    '-num_threads', threads,
    "-outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore'",
    '-out', file.path(out, 'host.txt')
)
system2('blastn', args = blastn.args)

# Read in the result ####
columns <- c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
             'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
hits <- read.delim(file.path(out, 'host.txt'), header = F, col.names = columns)

# Check to make sure only one primer hit was found and that the alignment contains the start of the primer ####
if(nrow(hits) < 1){
    
    stop('Primer does not align with the current settings.')
    
} else if(is.null(nrow(hits)) == T){
    
    stop('Primer does not align with the current settings.')
    
} else if(nrow(hits) > 1){
    
    stop('Primer aligns to multiple sites on the host sequence.')
    
} else if(hits$qstart != 1){
    
    stop('The start of the primer does not align to the host sequence.')
    
}

# If the host sequence is reversed, it needs to be flipped later ####
rev.seq <- F
if(hits$sstart > hits$send){
    
    rev.seq <- T
    
}

# Find the full length of the primer that will be blocked ####
primer.length <- readLines(list.files(in.path, pattern = '*primer.fa', full.names = T), n = 2)[2] |> nchar()

# Read in the host sequence as a vector ####
host.seq <- readLines(list.files(in.path, pattern = '*host.fa', full.names = T), n = 2)[2]
host.vect <- strsplit(host.seq, '')[[1]]

# Trim the vector and create the reverse complement if needed ####
if(rev.seq == T){
    
    host.trim <- host.vect |> head(hits$sstart)
    host.trim <- host.trim |> rev() |> paste(collapse = '') |>
        chartr('ATGCRYMKHDBV','TACGYRKMDHVB', x = _)
    host.trim <- strsplit(host.trim, '')[[1]]
    
} else {
    
    host.trim <- host.vect |> tail(-(hits$sstart) - 1)# |> paste(collapse = '')
    
}

# Define the minimum and maximum overlap of the blocking oligos ####
min.overlap <- 1
max.overlap <- ceiling(0.5 * primer.length)

# Determine the starting position and the number of shifts to investigate ####
starting <- primer.length - max.overlap + 1
oligo.length <- 30
oligo.indices <- starting:(starting + oligo.length - 1)
total.shifts <- max.overlap - min.overlap
index.shifts <- 0:total.shifts

# Obtain a table of candidate blocking oligos ####
block.out <- lapply(index.shifts, blocker.walker,
                    start = oligo.indices, vect = host.trim,
                    overlap = max.overlap) |> bind_rows()

# Write the sequences to a FASTA file ####
blockers <- block.out$seq
names(blockers) <- block.out$id
DNAStringSet(blockers) |> writeXStringSet(file.path(out, 'blockers.fa'))
