library(tidyverse)
library(qiime2R)


# Read in all of the necessary files
md <- read_tsv("sample-metadata.tsv")
ar_table <- read_qza("AR-rarefied-filtered-table.qza")
bray <- read_qza("AR-core-metrics-results/bray_curtis_pcoa_results.qza")
uUnifrac <- read_qza("AR-core-metrics-results/unweighted_unifrac_pcoa_results.qza")
bray_biofilm <- read_qza("biofilm-core-metrics-results/bray_curtis_pcoa_results.qza")
uUnifrac_biofilm <- read_qza("biofilm-core-metrics-results/unweighted_unifrac_pcoa_results.qza")

# bray curtis plot - ARs
bray_plt <- bray$data$Vectors %>%
  left_join(md) %>%
  ggplot(aes(x=PC1, y=PC2, color=Sample_Identifier)) +
  geom_point(size = 5) +
  xlab(paste("PC1: ", round(100*bray$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*bray$data$ProportionExplained[2]), "%")) +
  stat_ellipse(aes(PC1, PC2, group=Pipe_Material, color=Pipe_Material), type = "norm", linetype = 2) +
  theme_bw() +
  theme(legend.box = "horizontal", legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  ggsave("images/bray-PCA.png", width = 9.67, height = 4.46, units = "in", dpi = 300)


## Bray with Biofilm
bray_biofilm_plt <- bray_biofilm$data$Vectors %>%
  left_join(md) %>%
  ggplot(aes(x=PC1, y=PC2, color = Pipe_Material)) +
  geom_point(size = 5) +
  xlab(paste("PC1: ", round(100*bray$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*bray$data$ProportionExplained[2]), "%")) +
  stat_ellipse(aes(PC1, PC2, group=Pipe_Material, color=Pipe_Material), type = "norm", linetype = 2) +
  theme_bw() +
  theme(legend.box = "horizontal", legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  ggsave("images/bray-biofilm-PCA.png", width = 9.67, height = 4.46, units = "in", dpi = 300)

# uUnifrac plot
uUnifrac_plt <- uUnifrac$data$Vectors %>%
  left_join(md) %>%
  ggplot(aes(x=PC1, y=PC2, color=Sample_Identifier)) +
  geom_point(size = 5) +
  xlab(paste("PC1: ", round(100*uUnifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uUnifrac$data$ProportionExplained[2]), "%")) +
  stat_ellipse(aes(PC1, PC2, group=Pipe_Material, color=Pipe_Material), type = "norm", linetype = 2) +
  theme_bw() +
  theme(legend.box = "horizontal", legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  ggsave("images/uUnifrac-PCA.png", width = 9.67, height = 4.46, units = "in", dpi = 300)


## uUnifrac with Biofilm
uUnifrac_biofilm_plt <- uUnifrac_biofilm$data$Vectors %>%
  left_join(md) %>%
  ggplot(aes(x=PC1, y=PC2, color = Pipe_Material)) +
  geom_point(size = 5) +
  xlab(paste("PC1: ", round(100*uUnifrac_biofilm$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uUnifrac_biofilm$data$ProportionExplained[2]), "%")) +
  stat_ellipse(aes(PC1, PC2, group=Pipe_Material, color=Pipe_Material), type = "norm", linetype = 2) +
  theme_bw() +
  ggtitle("uUnifrac Distance PCA shows Pipe Biofilm are compositionally distinct from ARs")


# plot first distances ----------------------------------------------------

## uUnifrac Distance
# remove first row that gives column types
first_dist_bl <- read_tsv("longitudinal/first-distance-baseline-month2.tsv")[-1, ] %>% 
  select(id, Distance) %>% 
  left_join(md, by = c("id" = "SampleID")) %>% 
  mutate(Distance = as.numeric(Distance), 
         Start_Cal = factor(Start_Cal, levels = c("3-Sep", "4-Oct", "5-Nov", "6-Dec", "7-Jan", 
                                                  "8-Feb", "9-Mar", "10-Apr", "12-Jun", "13-Jul", 
                                                  "14-Aug", "15-Sep", "16-Oct")))

# generate boxplot
uUnifrac_dist_plt <- ggplot(first_dist_bl, aes(Start_Cal, Distance)) +
  geom_boxplot(outlier.color = "red") +
  facet_wrap(~ Pipe_Material) +
  scale_y_continuous(limits = c(0, 1.0)) + 
  ylab("uUnifrac Distance from First Month") +
  xlab("Months Since Experiment Start") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("images/uUnifrac-dist-baseline.png", width = 9.67, height = 4.46, units = "in", dpi = 300)



# plot distance between successive time points ----------------------------

first_dist <- read_tsv("longitudinal/first-distances.tsv")[-1, ] %>% 
  select(id, Distance) %>% 
  left_join(md, by = c("id" = "SampleID")) %>% 
  mutate(Distance = as.numeric(Distance), 
         Start_Cal = factor(Start_Cal, levels = c("3-Sep", "4-Oct", "5-Nov", "6-Dec", "7-Jan", 
                                                  "8-Feb", "9-Mar", "10-Apr", "12-Jun", "13-Jul", 
                                                  "14-Aug", "15-Sep", "16-Oct")),
         Pipe_Material = factor(Pipe_Material))

uUnifrac_first_dist_plt <- ggplot(first_dist, aes(Start_Cal, Distance)) +
  geom_boxplot(outlier.color = "red") +
  facet_wrap(~ Pipe_Material) +
  scale_y_continuous(limits = c(0, 1.0)) + 
  ylab("uUnifrac Distance Between Successive Time Points") +
  xlab("Months Since Experiment Start") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("images/uUnifrac-dist-successive.png", width = 9.67, height = 4.46, units = "in", dpi = 300)


## test for significance
library(DescTools)

kruskal.test(Distance ~ Pipe_Material, data = first_dist)

## test pairwise within groups
cast_iron_first_dist <- first_dist %>% 
  filter(Pipe_Material == 'Cast Iron')

cement_first_dist <- first_dist %>% 
  filter(Pipe_Material == 'Cement')

DunnTest(Distance ~ Start_Cal, data = cast_iron_first_dist, method = "fdr", alternative = "two.sided")
DunnTest(Distance ~ Start_Cal, data = cement_first_dist, method = "fdr", alternative = "two.sided")


# plot seasonal diferences in each AR -------------------------------------

cast_iron_season <- read_tsv("cast-iron-season-PERMANOVA.tsv") %>% select(-X1)
cement_season <- read_tsv("cement-season-PERMANOVA.tsv") %>% select(-X1)
