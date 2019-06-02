
# Plot Alpha Diversity ----------------------------------------------------

library(tidyverse)
library(qiime2R)
library(ggthemes)

# Read in all of the necessary files
md <- read_tsv("sample-metadata.tsv")
ar_table <- read_qza("AR-rarefied-table.qza")
taxonomy <- read_qza("taxonomy.qza")
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
  rename("Faith PD" = "faith_pd", "Observed OTUs" = "observed_otus",
           "Shannon" = "shannon", "Simpson" = "simpson")

# format long
df2 <- df %>% 
  gather(key = alpha_metric, value = value, `Shannon`:`Observed OTUs`)

# create plots ------------------------------------------------------------

alpha_plt <- ggplot(df2, aes(Start_Cal, value, color = Pipe_Material)) +
  geom_point(alpha = 0.2) +
  geom_smooth(aes(group = Pipe_Material)) +
  facet_wrap(~ alpha_metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  guides(colour = guide_legend(nrow = 1, title = "Pipe Material")) +
  ggtitle("Alpha Diversity Increases Over Time") +
  ggsave("images/all-alpha-diversity.png", width = 9.67, height = 4.46, units = "in", dpi = 600)
  

# alpha diversity by pipe material
p <- ggplot(df2, aes(x = Pipe_Material, y = value)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.color = "red") +
  facet_wrap(~alpha_metric, scales = "free_y") +
  theme_bw()
  

# alpha diversity by AR
p_samp_id <- ggplot(df2, aes(x = Sample_Identifier, y = value)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.color = "red") +
  facet_wrap(~alpha_metric, scales = "free_y") +
  theme_bw()

# Bonus Pairwise Testing --------------------------------------------------

df3 <- df %>% 
  mutate(Pipe_Material = factor(Pipe_Material, levels = unique(df$Pipe_Material)),
         Sample_Identifier = factor(Sample_Identifier, levels = unique(df$Sample_Identifier)))

kruskal.test(shannon ~ Pipe_Material, data = df3)
kruskal.test(simpson ~ Pipe_Material, data = df3)
kruskal.test(faith_pd ~ Pipe_Material, data = df3)
kruskal.test(shannon ~ Sample_Identifier, data = df3)
kruskal.test(simpson ~ Sample_Identifier, data = df3)
kruskal.test(faith_pd ~ Sample_Identifier, data = df3)

## pairwise tests - Dunn
library(DescTools)

DunnTest(shannon ~ Pipe_Material, data = df3, method = "fdr", alternative = "two.sided")
DunnTest(simpson ~ Pipe_Material, data = df3, method = "fdr", alternative = "two.sided")
DunnTest(faith_pd ~ Pipe_Material, data = df3, method = "fdr", alternative = "two.sided")
DunnTest(shannon ~ Sample_Identifier, data = df3, method = "fdr", alternative = "two.sided")
DunnTest(simpson ~ Sample_Identifier, data = df3, method = "fdr", alternative = "two.sided")
DunnTest(faith_pd ~ Sample_Identifier, data = df3, method = "fdr", alternative = "two.sided")


# plot alpha diversity over time ------------------------------------------

## Plot Faith's PD over Time
faith_months <- df2 %>% 
  filter(alpha_metric == "faith_pd") %>% 
  ggplot(aes(x = Start_Cal, y = value, color = Pipe_Material)) +
  stat_smooth(aes(group = Pipe_Material)) +
  geom_point(aes(color = Sample_Identifier), alpha = 0.2) +
  geom_line(aes(group = Sample_Identifier, color = Sample_Identifier), alpha = 0.2) +
  xlab("Months Since Experiment Start") +
  ylab("Faith's PD") +
  labs(color = "Sample Identifier") +
  theme_hc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.box = "horizontal") +
  guides(colour = guide_legend(nrow = 1)) +
  ggsave("images/faith-months.png", width = 9.67, height = 4.46, units = "in", dpi = 300)

## Shannon Index over time
shannon_months <- df2 %>% 
  filter(alpha_metric == "shannon") %>% 
  ggplot(aes(x = Start_Cal, y = value, color = Pipe_Material)) +
  stat_smooth(aes(group = Pipe_Material)) +
  geom_point(aes(color = Sample_Identifier), alpha = 0.2) +
  geom_line(aes(group = Sample_Identifier, color = Sample_Identifier), alpha = 0.2) +
  xlab("Months Since Experiment Start") +
  ylab("Shannon Index") +
  labs(color = "Sample Identifier") +
  theme_hc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.box = "horizontal") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(limits = c(2, 7)) +
  ggsave("images/shannon-months.png", width = 9.67, height = 4.46, units = "in", dpi = 300)


## Simpson Index over time
simpson_months <- df2 %>% 
  filter(alpha_metric == "simpson") %>% 
  ggplot(aes(x = Start_Cal, y = value, color = Pipe_Material)) +
  stat_smooth(aes(group = Pipe_Material)) +
  geom_point(aes(color = Sample_Identifier), alpha = 0.2) +
  geom_line(aes(group = Sample_Identifier, color = Sample_Identifier), alpha = 0.2) +
  xlab("Months Since Experiment Start") +
  ylab("Simpson Index") +
  labs(color = "Sample Identifier") +
  theme_hc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.box = "horizontal") +
  scale_y_continuous(limits = c(0.6, 1.0)) +
  guides(colour = guide_legend(nrow = 1)) +
  ggsave("images/simpson-months.png", width = 9.67, height = 4.46, units = "in", dpi = 300)

  
# track rate of change over time ------------------------------------------

faith_rate <- read_tsv("longitudinal/faith-first-differences.tsv")[-1, ]
shannon_rate <- read_tsv("longitudinal/shannon-first-differences.tsv")[-1, ]
simpson_rate <- read_tsv("longitudinal/simpson-first-differences.tsv") [-1, ]

## plot rate of change in Faith PD over months
faith_rate %>% 
  mutate(Months_Since_Start = as.numeric(Months_Since_Start),
         Difference = as.numeric(Difference)) %>% 
  ggplot(aes(x = Months_Since_Start, y = Difference, color=Sample_Identifier)) +
  geom_line(alpha = 0.2) +
  geom_smooth(aes(color=factor(Pipe_Material))) +
  geom_abline(slope = 0, ymin = 0, ymax = 0, color = "red", linetype = 2) +
  scale_x_continuous(breaks = 3:13) +
  scale_y_continuous(breaks = -2.5:5) +
  ggtitle("Rate of Change in Faith PD over time") +
  xlab("Months Since Experiment Start") +
  ylab("Difference in Faith PD from Previous Month") +
  theme_bw()


shannon_rate %>% 
  mutate(Months_Since_Start = as.numeric(Months_Since_Start),
         Difference = as.numeric(Difference)) %>% 
  ggplot(aes(x = Months_Since_Start, y = Difference, color=Sample_Identifier)) +
  geom_line(alpha = 0.2) +
  geom_smooth(aes(color=factor(Pipe_Material))) +
  geom_abline(slope = 0, ymin = 0, ymax = 0, color = "red", linetype = 2) +
  scale_x_continuous(breaks = 3:13) +
  scale_y_continuous(breaks = -2.5:5) +
  ggtitle("Rate of Change in Shannon Index over time") +
  xlab("Months Since Experiment Start") +
  ylab("Difference in Shannon Index from Previous Month") +
  theme_bw()

