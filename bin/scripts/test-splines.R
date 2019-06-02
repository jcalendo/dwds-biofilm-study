
# SplinectomeR analysis ---------------------------------------------------

library(qiime2R)
library(tidyverse)
library(splinectomeR)

# Read in all of the necessary files
md <- read_tsv("sample-metadata.tsv")
shannon <- read_qza("AR-core-metrics-results/shannon_vector.qza")
simpson <- read_qza("AR-core-metrics-results/simpson_vector.qza")
faith_pd <- read_qza("AR-core-metrics-results/faith_pd_vector.qza")
observed_otu <- read_qza("AR-core-metrics-results/observed_otus_vector.qza")

# Coerce alpha diversity metrics to a ggplot friendly data.frame
shannon_df <- shannon$data %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

simpson_df <- simpson$data %>% 
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

faith_df <- faith_pd$data %>% 
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

otu_df <-observed_otu$data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SampleID")

# join alpha metrics to metadata ------------------------------------------

df <- md %>% 
  left_join(shannon_df, by = "SampleID") %>% 
  left_join(simpson_df, by = "SampleID") %>% 
  left_join(faith_df, by = "SampleID") %>% 
  left_join(otu_df, by = "SampleID") %>% 
  mutate(Start_Cal = factor(Start_Cal, levels = c("2-Aug", "3-Sep", "4-Oct", "5-Nov", "6-Dec", "7-Jan", 
                                                  "8-Feb", "9-Mar", "10-Apr", "12-Jun", "13-Jul", 
                                                  "14-Aug", "15-Sep", "16-Oct"))) %>% 
  filter(Pipe_Material %in% c("Cast Iron", "Cement")) %>% 
  as.data.frame()


# permuspliner ------------------------------------------------------------

# check to see if there is a greater-than-chance between two groups of interest
# over the x-variable

faith_result <- permuspliner(data = df, x = 'Days_Since_Start', y = 'faith_pd',
                       cases = 'Sample_Identifier', category = 'Pipe_Material', 
                       groups = c('Cast Iron','Cement'),  cut_sparse = 2, retain_perm = TRUE)

## test Shannon Diversity
shannon_result <- permuspliner(data = df, x = 'Days_Since_Start', y = 'shannon',
                               cases = 'Sample_Identifier', category = 'Pipe_Material', 
                               groups = c('Cast Iron','Cement'),  cut_sparse = 2, retain_perm = TRUE)

# test simpson diversity
simpson_result <- permuspliner(data = df, x = 'Days_Since_Start', y = 'simpson',
                               cases = 'Sample_Identifier', category = 'Pipe_Material', 
                               groups = c('Cast Iron','Cement'),  cut_sparse = 2, retain_perm = TRUE)

otus_result <- permuspliner(data = df, x = 'Days_Since_Start', y = 'observed_otus',
                            cases = 'Sample_Identifier', category = 'Pipe_Material', 
                            groups = c('Cast Iron','Cement'),  cut_sparse = 2, retain_perm = TRUE)

# test for non-zero slope -------------------------------------------------
cast_iron_faith_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'faith_pd',
                        cases = 'Sample_Identifier', category = 'Pipe_Material', 
                        group = 'Cast Iron', perms = 999)

cement_faith_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'faith_pd',
                                       cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                    group = 'Cement', perms = 999)

cast_iron_shannon_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'shannon',
                                         cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                         group = 'Cast Iron', perms = 999)

cement_shannon_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'shannon',
                                         cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                         group = 'Cement', perms = 999)

cast_iron_simpson_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'simpson',
                                         cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                         group = 'Cast Iron', perms = 999)

cement_simpson_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'simpson',
                                      cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                      group = 'Cement', perms = 999)

cast_iron_otu_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'observed_otus',
                                     cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                     group = 'Cast Iron', perms = 999)

cement_otu_trend <- trendyspliner(data = df, x = 'Days_Since_Start', y = 'observed_otus',
                                     cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                     group = 'Cement', perms = 999)

# test for significance at any X ------------------------------------------

# test whether two groups differ at any point along the x-axis even though the overall 
# difference may not be significant
faith_sliding <- sliding_spliner(data = df, xvar = 'Days_Since_Start', yvar = 'shannon',
                          category = 'Pipe_Material', groups = c('Cast Iron','Cement'), 
                          cases = 'Sample_Identifier', ints = 15, cut_low = 2)

shannon_sliding <- sliding_spliner(data = df, xvar = 'Days_Since_Start', yvar = 'shannon',
                                 category = 'Pipe_Material', groups = c('Cast Iron','Cement'), 
                                 cases = 'Sample_Identifier', ints = 15, cut_low = 2)

simpson_sliding <- sliding_spliner(data = df, xvar = 'Days_Since_Start', yvar = 'simpson',
                                   category = 'Pipe_Material', groups = c('Cast Iron','Cement'), 
                                   cases = 'Sample_Identifier', ints = 15, cut_low = 2)

otus_sliding <- sliding_spliner(data = df, xvar = 'Days_Since_Start', yvar = 'observed_otus',
                                category = 'Pipe_Material', groups = c('Cast Iron','Cement'), 
                                cases = 'Sample_Identifier', ints = 15, cut_low = 2)

# create prettier plots ---------------------------------------------------

## make the faith_sliding results a little prettier
faith_sliding_plt <- faith_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 0.5)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("There are no significant differences in Faith PD between Pipe Materials at any time") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

shannon_sliding_plt <- shannon_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 0.5)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("There are no significant differences in Shannon Index between Pipe Materials at any time") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

simpson_sliding_plt <- simpson_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 0.8)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("There are no significant differences in Simpson Diversity between Pipe Materials at any time") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()


# test feature abundance splines ------------------------------------------

features <- read_tsv("family-level-feature-abundance.tsv")[-1, ]

## Mycobacteriaceae

myco <- features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Corynebacteriales;D_4__Mycobacteriaceae`) %>% 
  rename("Mycobacteriaceae" = "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Corynebacteriales;D_4__Mycobacteriaceae") %>%
  mutate(Mycobacteriaceae = as.numeric(Mycobacteriaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
myco_perm <- permuspliner(data = myco, x = 'Months_Since_Start', y = 'Mycobacteriaceae',
                          cases = 'Sample_Identifier', category = 'Pipe_Material', 
                          groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron
cast_iron_myco_trend <- trendyspliner(data = myco, x = 'Months_Since_Start', y = 'Mycobacteriaceae',
                                       cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                       group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_myco_trend <- trendyspliner(data = myco, x = 'Months_Since_Start', y = 'Mycobacteriaceae',
                                      cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                      group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
myco_sliding <- sliding_spliner(data = myco, xvar = "Months_Since_Start", yvar = "Mycobacteriaceae",
                                 category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                 cases = "Sample_Identifier", cut_low = 1, ints = 15)

myco_sliding_plt <- myco_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 0.8)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Mycobacteriaceae Relative Abundances") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Legionellaceae

lego <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Legionellales;D_4__Legionellaceae`) %>% 
  rename("Legionellaceae" = "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Legionellales;D_4__Legionellaceae") %>%
  mutate(Legionellaceae = as.numeric(Legionellaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
lego_perm <- permuspliner(data = lego, x = 'Months_Since_Start', y = 'Legionellaceae',
                          cases = 'Sample_Identifier', category = 'Pipe_Material', 
                          groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_lego_trend <- trendyspliner(data = lego, x = 'Months_Since_Start', y = 'Legionellaceae',
                                      cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                      group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_lego_trend <- trendyspliner(data = lego, x = 'Months_Since_Start', y = 'Legionellaceae',
                                   cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                   group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
lego_sliding <- sliding_spliner(data = lego, xvar = "Months_Since_Start", yvar = "Legionellaceae",
                                category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
lego_sliding_plt <- lego_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 0.8)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Legionellaceae Relative Abundances") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Burkholderiaceae

burk <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Burkholderiaceae`) %>% 
  rename("Burkholderiaceae" = "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Burkholderiaceae") %>%
  mutate(Burkholderiaceae = as.numeric(Burkholderiaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
burk_perm <- permuspliner(data = burk, x = 'Months_Since_Start', y = 'Burkholderiaceae',
                          cases = 'Sample_Identifier', category = 'Pipe_Material', 
                          groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_burk_trend <- trendyspliner(data = burk, x = 'Months_Since_Start', y = 'Burkholderiaceae',
                                      cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                      group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_burk_trend <- trendyspliner(data = burk, x = 'Months_Since_Start', y = 'Burkholderiaceae',
                                   cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                   group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
burk_sliding <- sliding_spliner(data = burk, xvar = "Months_Since_Start", yvar = "Burkholderiaceae",
                                category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
burk_sliding_plt <- burk_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Burkholderiaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Azospirillaceae

azo <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Azospirillales;D_4__Azospirillaceae`) %>% 
  rename("Azospirillaceae" = "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Azospirillales;D_4__Azospirillaceae") %>%
  mutate(Azospirillaceae = as.numeric(Azospirillaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
azo_perm <- permuspliner(data = azo, x = 'Months_Since_Start', y = 'Azospirillaceae',
                          cases = 'Sample_Identifier', category = 'Pipe_Material', 
                          groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_azo_trend <- trendyspliner(data = azo, x = 'Months_Since_Start', y = 'Azospirillaceae',
                                      cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                      group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_azo_trend <- trendyspliner(data = azo, x = 'Months_Since_Start', y = 'Azospirillaceae',
                                   cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                   group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
azo_sliding <- sliding_spliner(data = azo, xvar = "Months_Since_Start", yvar = "Azospirillaceae",
                                category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
azo_sliding_plt <- azo_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Azospirillaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Nitrosomonadaceae

nitrosomo <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Nitrosomonadaceae`) %>% 
  rename("Nitrosomonadaceae" = "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Nitrosomonadaceae") %>%
  mutate(Nitrosomonadaceae = as.numeric(Nitrosomonadaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
nitrosomo_perm <- permuspliner(data = nitrosomo, x = 'Months_Since_Start', y = 'Nitrosomonadaceae',
                         cases = 'Sample_Identifier', category = 'Pipe_Material', 
                         groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_nitrosomo_trend <- trendyspliner(data = nitrosomo, x = 'Months_Since_Start', y = 'Nitrosomonadaceae',
                                     cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                     group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_nitrosomo_trend <- trendyspliner(data = nitrosomo, x = 'Months_Since_Start', y = 'Nitrosomonadaceae',
                                  cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                  group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
nitrosomo_sliding <- sliding_spliner(data = nitrosomo, xvar = "Months_Since_Start", yvar = "Nitrosomonadaceae",
                               category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                               cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
nitrosomo_sliding_plt <- nitrosomo_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Nitrosomonadaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Nitrospiraceae

nitrospira <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Nitrospirae;D_2__Nitrospira;D_3__Nitrospirales;D_4__Nitrospiraceae`) %>% 
  rename("Nitrospiraceae" = "D_0__Bacteria;D_1__Nitrospirae;D_2__Nitrospira;D_3__Nitrospirales;D_4__Nitrospiraceae") %>%
  mutate(Nitrospiraceae = as.numeric(Nitrospiraceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
nitrospira_perm <- permuspliner(data = nitrospira, x = 'Months_Since_Start', y = 'Nitrospiraceae',
                               cases = 'Sample_Identifier', category = 'Pipe_Material', 
                               groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_nitrospira_trend <- trendyspliner(data = nitrospira, x = 'Months_Since_Start', y = 'Nitrospiraceae',
                                           cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                           group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_nitrospira_trend <- trendyspliner(data = nitrospira, x = 'Months_Since_Start', y = 'Nitrospiraceae',
                                        cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                        group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
nitrospira_sliding <- sliding_spliner(data = nitrospira, xvar = "Months_Since_Start", yvar = "Nitrospiraceae",
                                     category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                     cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
nitrospira_sliding_plt <- nitrospira_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Nitrospiraceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Pseudomonadaceae

pseudo <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Nitrospirae;D_2__Nitrospira;D_3__Nitrospirales;D_4__Nitrospiraceae`) %>% 
  rename("Pseudomonadaceae" = "D_0__Bacteria;D_1__Nitrospirae;D_2__Nitrospira;D_3__Nitrospirales;D_4__Nitrospiraceae") %>%
  mutate(Pseudomonadaceae = as.numeric(Pseudomonadaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
pseudo_perm <- permuspliner(data = pseudo, x = 'Months_Since_Start', y = 'Pseudomonadaceae',
                                cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_pseudo_trend <- trendyspliner(data = pseudo, x = 'Months_Since_Start', y = 'Pseudomonadaceae',
                                            cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                            group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_pseudo_trend <- trendyspliner(data = pseudo, x = 'Months_Since_Start', y = 'Pseudomonadaceae',
                                         cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                         group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
pseudo_sliding <- sliding_spliner(data = pseudo, xvar = "Months_Since_Start", yvar = "Pseudomonadaceae",
                                      category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                      cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
pseudo_sliding_plt <- pseudo_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Pseudomonadaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Rhodocyclaceae

rhodo <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Rhodocyclaceae`) %>% 
  rename("Rhodocyclaceae" = "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Rhodocyclaceae") %>%
  mutate(Rhodocyclaceae = as.numeric(Rhodocyclaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
rhodo_perm <- permuspliner(data = rhodo, x = 'Months_Since_Start', y = 'Rhodocyclaceae',
                            cases = 'Sample_Identifier', category = 'Pipe_Material', 
                            groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_rhodo_trend <- trendyspliner(data = rhodo, x = 'Months_Since_Start', y = 'Rhodocyclaceae',
                                        cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                        group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_rhodo_trend <- trendyspliner(data = rhodo, x = 'Months_Since_Start', y = 'Rhodocyclaceae',
                                     cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                     group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
rhodo_sliding <- sliding_spliner(data = rhodo, xvar = "Months_Since_Start", yvar = "Rhodocyclaceae",
                                  category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                  cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
rhodo_sliding_plt <- rhodo_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Rhodocyclaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()

## Chitinophagaceae

chit <-features %>% 
  select(Sample_Identifier, Pipe_Material, Months_Since_Start, `D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae`) %>% 
  rename("Chitinophagaceae" = "D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae") %>%
  mutate(Chitinophagaceae = as.numeric(Chitinophagaceae),
         Months_Since_Start = as.numeric(Months_Since_Start)) %>% 
  as.data.frame()

# test overall significance
chit_perm <- permuspliner(data = chit, x = 'Months_Since_Start', y = 'Chitinophagaceae',
                           cases = 'Sample_Identifier', category = 'Pipe_Material', 
                           groups = c('Cast Iron','Cement'),  cut_sparse = 1, retain_perm = TRUE)

# test for significance departure from null slope for Cast Iron - p = 0.001
cast_iron_chit_trend <- trendyspliner(data = chit, x = 'Months_Since_Start', y = 'Chitinophagaceae',
                                       cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                       group = 'Cast Iron', perms = 999)

# test for significant departure from null slope for Cement ARs - p = 0.002
cement_chit_trend <- trendyspliner(data = chit, x = 'Months_Since_Start', y = 'Chitinophagaceae',
                                    cases = 'Sample_Identifier', category = 'Pipe_Material', 
                                    group = 'Cement', perms = 999)

# test for significant difference in abundance at any X
chit_sliding <- sliding_spliner(data = chit, xvar = "Months_Since_Start", yvar = "Chitinophagaceae",
                                 category = "Pipe_Material", groups = c('Cast Iron','Cement'), 
                                 cases = "Sample_Identifier", cut_low = 1, ints = 15)

# create pvalue plot
chit_sliding_plt <- chit_sliding$pval_table %>% 
  ggplot(aes(Months_Since_Start, p_value)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.05, color = "red", linetype = 2) +
  scale_y_continuous(limits = c(0.01, 1.0)) +
  scale_x_continuous(breaks = seq(2, 12, 1)) +
  ggtitle("No Significant Differences in Chitinophagaceae Relative Abundance Curves") +
  xlab("Months Since Experiment Start") +
  ylab("p_value") +
  theme_bw()