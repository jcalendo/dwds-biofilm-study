
# calculate and plot correlations -----------------------------------------

library(tidyverse)
library(Hmisc)
library(corrplot)


# read in dataframe with alpha diversity metrics and
# eliminate non-numeric columns 
alpha_metrics <- read_tsv("alpha_diversity_metrics.tsv") %>% 
  select(SampleID, Pipe_Material, PWD_Inf_Total_Chlorine:Days_Since_Start, Avg_Temp_Inf,
         Months_Since_Start, Months_Since_Start_Exact, alpha_metric,
         value) %>% 
  spread(key = alpha_metric, value = `value`) %>% 
  column_to_rownames(var = 'SampleID')


# clean up column names for plotting
alpha_metrics2 <- alpha_metrics %>% 
  rename("Total Chlorine" = "PWD_Inf_Total_Chlorine",
         "Total Conductivity" = "PWD_Inf_Total_Conductivity",
         "Temperature" = "PWD_Inf_Temp",
         "Turbidity" = "PWD_Inf_Turb",
         "Max Chlorine" = "Max_PWD_Inf_Total_Chlorine",
         "Max Conductivity" = "Max_PWD_Inf_Total_Conductivity",
         "Max Temperature" = "Max_PWD_Inf_Temp",
         "Max Turbidity" = "Max_PWD_Inf_Turb",
         "Months Since Start" = "Months_Since_Start",
         "Faith PD" = "faith_pd",
         "Shannon" = "shannon",
         "Simpson" = "simpson") %>% 
  select(-Days_Since_Start, -Months_Since_Start_Exact, -Avg_Temp_Inf)

# split by Pipe_material
cast_iron <- alpha_metrics2 %>% filter(Pipe_Material == 'Cast Iron') %>% select(-Pipe_Material)
cement <- alpha_metrics2 %>% filter(Pipe_Material == 'Cement') %>% select(-Pipe_Material)

# compute correlation -----------------------------------------------------

cast_iron_res <- rcorr(x = as.matrix(cast_iron), type = 'spearman')
cement_res <- rcorr(x = as.matrix(cement), type = 'spearman')

# plot results in correlation plot
corrplot(cast_iron_res$r, type="upper", 
         p.mat = cast_iron_res$P, sig.level = 0.05, insig = "blank")

corrplot(cement_res$r, type="upper",
                       p.mat = cement_res$P, sig.level = 0.05, insig = "blank")


