
# create krona plots from taxa barplots output ----------------------------

library(tidyverse)


taxa_levels <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
taxa <- read_csv("~/Downloads/level-7.csv") %>% 
  gather(key = Taxon, value = count, "D_0__Bacteria;D_1__Acidobacteria;D_2__Acidobacteriia;D_3__Solibacterales;D_4__Solibacteraceae (Subgroup 3);D_5__Bryobacter;D_6__metagenome":"D_0__Bacteria;D_1__WPS-2;D_2__metagenome;D_3__metagenome;D_4__metagenome;D_5__metagenome;D_6__metagenome") %>%
  separate(Taxon, into = taxa_levels, sep = ';') %>% 
  mutate_at(vars(taxa_levels), funs(str_replace(., "D_.__", ""))) %>% 
  mutate_at(vars(taxa_levels), funs(str_replace(., "__", ""))) %>% 
  filter(count > 0)

# create krona output format - cast iron
cast_iron <- taxa %>% 
  filter(Pipe_Material == 'Cast Iron') %>% 
  select(count, taxa_levels)

# create krona output for cement
cement <- taxa %>% 
  filter(Pipe_Material == 'Cement') %>% 
  select(count, taxa_levels)


write_tsv(cast_iron, "qiime2Krona/cast-iron.tsv")
write_tsv(cement, "qiime2Krona/cement.tsv")
