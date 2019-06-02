# Clean Metadata

library(tidyverse)
library(lubridate)

md <- read_tsv("sample-metadata.tsv")

# DEFINE START_DATE
# know from Dr. Murphy's meatdata file that 8/8/2017 was 49 days from the start
day_49 <- ymd("2017-08-08")
START_DATE = day_49 - days(49)


md2 <-md %>% 
  mutate(Days_Since_Start = time_length(START_DATE %--% Date, unit = "days"),
         Months_Since_Start_Exact = time_length(START_DATE %--% Date, unit = "month"),
         Calendar_Month = month(Date, label = TRUE),
         Start_Cal = str_c(Months_Since_Start, Calendar_Month, sep = "-"))


write_tsv(md2, "sample-metadata.tsv")