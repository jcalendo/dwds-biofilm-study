
# plot selected family level taxa -----------------------------------------
library(tidyverse)


level_5 <- read_csv("family-level-differentially-abundant-taxa.csv")

# gather level_5 table into long form
taxa_df <- level_5 %>% 
  gather(key = taxa, value = count, "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Kineosporiales;D_4__Kineosporiaceae":"D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Rhodocyclaceae") %>% 
  mutate(family = str_extract(taxa, "D_4__.*")) %>% 
  mutate(family = str_replace(family, "D_4__", ""))

# create a filled plot
ggplot(taxa_df, aes(Date, count, fill = family)) +
  geom_area(position = 'fill') +
  facet_wrap(~ Pipe_Material) +
  theme_bw() +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~ Sample_Identifier) +
  xlab(NULL) +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("images/differentially-abundant-taxa.png", width = 9.67, height = 4.46, units = "in", dpi = 300)


# plot relative frequencies at family level -------------------------------

rel_freq_family <- read_tsv("family-level-relative-frequency-table/feature-table.tsv")
important_features <- read_tsv("ANCOM-features-relative-frequencies/feature-table.tsv") %>% distinct(FeatureID) %>% pull()
md <- read_tsv("sample-metadata.tsv")

# prepare for plotting
taxa_df2 <- rel_freq_family %>% 
  gather(key = SampleID, value = relative_abundance, `AR_1_1-19-18`:`AR-6-9_20_18`) %>% 
  filter(FeatureID %in% important_features) %>% 
  mutate(family = str_extract(FeatureID, "D_4__.*")) %>% 
  mutate(family = str_replace(family, "D_4__", "")) %>% 
  select(-FeatureID) %>% 
  left_join(md, by = "SampleID") 
  
# create plot
ggplot(taxa_df2, aes(Date, relative_abundance, fill = family)) +
  geom_area() +
  ylab('relative abd.') +
  xlab(NULL) +
  scale_y_continuous(limits = c(0, 0.5)) +
  facet_wrap(~ Pipe_Material, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  scale_fill_brewer(palette = "Spectral") +
  ggsave("images/relative-abd-family-taxa.png", width = 9.67, height = 4.46, units = "in", dpi = 300)



