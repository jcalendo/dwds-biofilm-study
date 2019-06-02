
# explore families --------------------------------------------------------

library(tidyverse)


taxa_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa <- read_csv("~/Downloads/level-7.csv") %>% 
  gather(key = taxon, 
         value = count, 
         "D_0__Bacteria;D_1__Acidobacteria;D_2__Acidobacteriia;D_3__Solibacterales;D_4__Solibacteraceae (Subgroup 3);D_5__Bryobacter;D_6__metagenome":"D_0__Bacteria;D_1__WPS-2;D_2__metagenome;D_3__metagenome;D_4__metagenome;D_5__metagenome;D_6__metagenome") %>% 
  separate(taxon, into = taxa_levels, sep = ";") %>% 
  mutate_at(vars(taxa_levels), funs(str_replace(., "D_[0-9]__", "")))

# define a function for estracting unique genera given family name
get_genera = function(fam) {
  result <- taxa %>% 
    filter(family == fam) %>% 
    group_by(genus) %>% 
    count()
  
  result
}

burk <- get_genera("Burkholderiaceae")
hypho <- get_genera(("Hyphomicrobiaceae"))
rhodo <- get_genera("Rhodocyclaceae")
sphingo <- get_genera(("Sphingomonadaceae"))


# find some pathogens! ----------------------------------------------------

# create vector of genera that contain common pathogens
pathogenic_genera <- c('Bacillus', 'Bartonella', 'Bordetella', 'Borrelia',
                       'Brucella', 'Campylobacter', 'Chlamydia', 'Chlamydophila',
                       'Clostridium', 'Corynebacterium', 'Enterococcus', 'Escherichia',
                       'Francisella', 'Haemophilus', 'Helicobacter', 'Legionella',
                       'Leptospira', 'Listeria', 'Mycobacterium', 'Mycoplasma', 
                       'Neisseria', 'Pseudomonas', 'Rickettsia', 'Salmonella', 
                       'Shigella', 'Staphylococcus', 'Streptococcus', 'Treponema',
                       'Ureaplasm', 'Vibrio', 'Yersinia')

pathogens <- taxa %>% 
  filter(genus %in% pathogenic_genera) %>% 
  filter(count > 0) %>% 
  group_by(genus, species) %>% 
  count()

psuedo <- taxa %>% 
  filter(genus == 'Pseudomonas' & count > 0) %>% 
  arrange(Date)


