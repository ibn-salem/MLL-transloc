# convert TSV file from JASPAR to BED file

require(tidyverse)

infile = "data/JASPAR2018/MA0139.1.tsv"

col_names = c("chr", "start", "end", "name", "score", "log10_pval_times_100", "strand")

df = read_tsv(infile, col_names = col_names, skip = 1, 
              col_type = cols(
                chr = col_character(),
                start = col_integer(),
                end = col_integer(),
                name = col_character(),
                score = col_integer(),
                log10_pval_times_100 = col_integer(),
                strand = col_character()
              ))

bed <- df %>% 
  # covnert from 1-based to 0-based half open
  mutate(start = start - 1) %>% 
  select(chr, start, end, name, log10_pval_times_100, strand)

write_tsv(bed, str_c(infile, ".bed"), col_names = FALSE)
