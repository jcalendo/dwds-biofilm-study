# Create Krona Plots in R

## load tidyverse and qiime2R libraries and read in artifacts and metadata
library(tidyverse)
library(qiime2R)


# load metadata
md <- read_tsv("sample-metadata.tsv")

# load feature table
ft <- read_qza("decontam-taxa-filtered-table.qza")

# load taxonomy
tax <- read_qza("taxonomy.qza")

## Extract the data from .qza artifacts and coerce to data frames
# extract data from feature table
ft_df <- ft$data %>% as.data.frame() %>% rownames_to_column(var = "Feature.ID")

# extract data from taxonomy file - change factor types to characters
tax_df <- tax$data %>% as_tibble(rownames = NULL) %>% mutate_if(is.factor, as.character)


## Join the taxonomy annotations to the feature table
annotated_ft <- ft_df %>% 
  left_join(tax_df, by = "Feature.ID") %>% 
  select(-Confidence)

## Reshape the table so it can be joined to metadata

# Once the table is in its long form the rows that contain counts of zero can simply be removed because they do not contain any relevant information for the krona plot.
long_annotated_ft <- annotated_ft %>% 
  gather(key = "SampleID", value = "count", `AR_1_1-19-18`:`MEB-11_16_18`) %>% 
  select(SampleID, count, Taxon, -Feature.ID) %>% 
  filter(count > 0)

## Join the metadata to the long table
complete_df <- long_annotated_ft %>% 
  left_join(md, by = c("SampleID" = "SampleID"))

## Split up the Taxon annotations into separate columns
# define a taxonomy label vector
taxa_labels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# split the Taxon column into separate columns and clean up the taxon labels
split_taxon_df <- complete_df %>% 
  separate(Taxon, into = taxa_labels, sep = ";") %>% 
  mutate_at(vars(taxa_labels), str_replace, pattern = ".*.__", replacement = "") %>% 
  replace(is.na(.), "")

## Filter on metadata to select samples of interest
### Select Cast Iron ARs only
# Once we have filtered the table on our metadata of interest we can drop all of the metadata columns that Krona does not need and write a tsv file that can be used as input to ktImportText
cast_iron_AR <- split_taxon_df %>% 
  filter(Pipe_Material == "Cast Iron" & Sample_Type == "AR") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/Cast-Iron-AR-krona.tsv", col_names = FALSE)

### Select only Cement ARs
cement_ARs <- split_taxon_df %>% 
  filter(Pipe_Material == "Cement" & Sample_Type == "AR") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/Cement-AR-krona.tsv", col_names = FALSE)

### Select only Pipe Biofilm
pipe_biofilm <- split_taxon_df %>% 
  filter(Pipe_Material == "Pipe Biofilm") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/Pipe-Biofilm-krona.tsv", col_names = FALSE)

### Bulk Water
bulk_water <- split_taxon_df %>% 
  filter(Pipe_Material == "Bulk-Water") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/Bulk-Water-krona.tsv", col_names = FALSE)

## Now create files for individual ARs
ar_1 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-1") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-1-krona.tsv", col_names = FALSE)

ar_2 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-2") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-2-krona.tsv", col_names = FALSE)

ar_3 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-3") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-3-krona.tsv", col_names = FALSE)

ar_4 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-4") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-4-krona.tsv", col_names = FALSE)

ar_5 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-5") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-5-krona.tsv", col_names = FALSE)

ar_6 <- split_taxon_df %>% 
  filter(Sample_Identifier == "AR-6") %>% 
  select(count, Kingdom:Species) %>% 
  write_tsv("qiime2Krona/AR-6-krona.tsv", col_names = FALSE)