
install.packages("tidyverse")

library(ggplot2)
library(tidyverse)
library(dplyr)

virus_data<-read.csv("/cloud/project/question-5-data/Cui_etal2014.csv")
  
virus_clean$Virion.volume..nm.nm.nm.<- as.numeric(virus_clean$Virion.volume..nm.nm.nm.)
  
colnames(virus_clean)[10] = "virion_volume"
colnames(virus_clean)[12] = "genome_length"

virus_clean <- mutate(virus_clean, log_virion_volume = log(virion_volume))
virus_clean <- mutate(virus_clean, log_genome_length = log(genome_length))

# this applies the log transformation to the data and adds these as new columns in the dataset


