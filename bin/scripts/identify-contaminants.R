# IDENTIFY CONTAMINANTS WITH DECONTAM
################################################################################
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(decontam)


# create phyloseq object from the qza artifacts
physeq <- qza_to_phyloseq(features = "table.qza", taxonomy = "taxonomy.qza", metadata = "sample-metadata.tsv")

# put sample_data into ggplot friendly df
df <- as.data.frame(sample_data(physeq)) 
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Plot Library Size for True and Control Samples
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

# Only include AR and MEB samples - Pipe Biofilm & DS are too variable
physeq <- subset_samples(physeq, Sample_Type %in% c("AR", "MEB"))
df <- as.data.frame(sample_data(physeq))
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Show new plot with only AR and MEB samples
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

# decontamination via prevalance
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

# list the possible contaminants by abundances
which(contamdf.prev$contaminant)

# identify using the taxon table
poss_contaminants <- tax_table(physeq)[which(contamdf.prev$contaminant)]

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

# plot presence-absence results                     
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# write the possible contaminants to a separate tsv file - used later to filter table
contams <- which(contamdf.prev$contaminant)
write.csv(tax_table(physeq)[contams], "contaminants.tsv")

contaminants <- read_csv("contaminants.csv") %>% write_tsv("contaminants.tsv")
