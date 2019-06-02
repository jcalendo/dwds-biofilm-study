library(tidyverse)
library(qiime2R)
library(phyloseq)
library(NMIT)
library(vegan)


# import data as phyloseq object ------------------------------------------

ar_table <- read_qza("AR-filtered-table.qza")
taxonomy <- read_qza("taxonomy.qza")
tree <- read_qza("rooted-tree.qza")
metadata <- read_tsv("sample-metadata.tsv")

# reshape the Taxonomy table
taxtable <- taxonomy$data %>% 
  as.tibble() %>% 
  mutate(Taxon = str_replace_all(Taxon, "D_.__", "")) %>% 
  separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

# create the phyloseq object
physeq <- phyloseq(
  otu_table(ar_table$data, taxa_are_rows = T),
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))



# perform PERMANOVA -------------------------------------------------------

physeq_bray <- phyloseq::distance(physeq, method = "bray")

# make datframe from sample_dta
sampledf <- data.frame(sample_data(physeq))

# adonis test
adonis(physeq_bray ~ (Pipe_Material / Sample_Identifier) + Start_Cal, data = sampledf, strata = sampledf$Pipe_Material)

# perform nmit ------------------------------------------------------------

nmit_result <- NMIT_phyloseq(physeq, id.var = "Sample_Identifier", cov.var = "Pipe_Material", time.var = "Months_Since_Start")

