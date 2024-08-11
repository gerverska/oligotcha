# Colorblind palette ####
blind <- c('#000000', '#E69F00', '#56B4E9',
           '#009E73', '#F0E442', '#0072B2',
           '#D55E00', '#CC79A7')

extended <- c('#fcf8ae',
              '#2bc6ff',
              '#ba9863',
              '#20a189',
              '#3c4278',
              '#cc1fb2',
              '#00ffbf',
              '#1a7878',
              '#3d2b0d',
              '#0025ba',
              '#ebb217',
              '#d4d4ff',
              '#027d3e',
              '#f2ff00',
              '#c642ff',
              '#ff0000',
              '#5252f7',
              '#ffa1a1',
              '#1e0078',
              '#780000',
              '#000000')

# ggplot2 theme ####
base.theme <- ggplot2::theme(text = ggplot2::element_text(size = 14),
                    axis.line = ggplot2::element_line(colour = 'black'),
                    axis.title.x = ggplot2::element_text(face = 'bold'),
                    axis.title.y = ggplot2::element_text(face = 'bold'),
                    legend.title = ggplot2::element_text(face = 'bold'))

# 01-design ####
# Compile blocking oligo candidates ####
blocker.walker <- function(i, start, vect, overlap){
    
    shifted.indices <- start + i
    shifted.seq <- vect[shifted.indices] |> paste(collapse = '')
    
    current.overlap <- max.overlap - i
    
    id <- paste0('base', shifted.indices[[1]], '_overlap', current.overlap)
    
    data.frame(id = id, seq = shifted.seq)
    
}

# 02-analyze.r ####
# Parse taxonomy ####
tax <- function(x){
    
    x$stitle <- x$stitle |> strsplit(' ') |> sapply( '[', 2)
    
    x$kingdom <- str_extract(x$stitle, "Fungi;") |> str_remove(';')
    x$phylum <- str_extract(x$stitle, "[[:alpha:]]+mycota;") |> str_remove(';')
    x$class <- str_extract(x$stitle, "[[:alpha:]]+mycetes;") |> str_remove(';')
    x$order <- str_extract(x$stitle, "[[:alpha:]]+ales;") |> str_remove(';')
    x$family <- str_extract(x$stitle, "[[:alpha:]]+aceae;") |> str_remove(';')
    
    x
    
}