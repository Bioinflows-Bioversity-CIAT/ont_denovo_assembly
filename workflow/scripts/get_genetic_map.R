rm(list = ls())

pks <- c('dplyr','argparse')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', names(pks)[!pks], ' is/are not available for load.')

pr <- ArgumentParser()

pr$add_argument("-map", type = "character", metavar = 'file',
                help = "Path to genetic map")

pr$add_argument("-blast", type = "character", metavar = "file", 
                help = "Path to blastp (fmt6) output file")

pr$add_argument("-out", type = "character", metavar = "path",
                help = "Outputh path")                    

ar = pr$parse_args()

# Functions

# Function to merge the genetic map and the blast
# output of the genetic markers and the genome draft
get_genetic_map <- function(map_path, blast_path){
  map <- read.table(map_path, sep = '\t', header = F)
  colnames(map) <- c("Chr", "marker_id", "gen_pos")
  
  blast <- read.table(blast_path, sep = '\t', header = F)
  colnames(blast) <-  c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
  
  blast <- blast %>% 
    group_by(qseqid) %>% 
    mutate(avg_pos = round((sstart + send)/2, 0)) %>% 
    rename(marker_id = qseqid)
  
  out <- merge(blast, map, by = "marker_id") %>% 
    select(sseqid, avg_pos, Chr, gen_pos) %>% 
    mutate(pos_e = avg_pos + 1,
      lg_tag = glue::glue("{Chr}:{gen_pos}")) %>%
    rename(
      contig_id = sseqid,
      marker_pos = avg_pos,
      lg = Chr) %>%
    arrange(lg, gen_pos) %>%
    select(contig_id, marker_pos, pos_e, lg_tag)
}

# Function call

out <- get_genetic_map(map_path = ar$map,
                blast_path = ar$blast)

write.table(out, ar$out, sep = '\t', row.names = F, col.names = F, quote = F)