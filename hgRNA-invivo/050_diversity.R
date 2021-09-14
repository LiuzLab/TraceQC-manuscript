library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

all_data <- readRDS("./data/040.2_all_data_normed.rds")

############## calculate measures ##############
info_table <- all_data$info_table
data_summary <- data.frame()
for (sample_name in names(all_data)[2:length(all_data)]) {
  tmp <- filter(info_table,Run==sample_name)

  sample <- all_data[[sample_name]] %>%
    group_by(identifier,type,start,length,mutate_to) %>%
    summarise(count=sum(count)) %>%
    ungroup

  shannon_entropy <- group_by(sample,identifier) %>%
    summarise(SE=Entropy(count))
    
  unmutated <- all_data[[sample_name]] %>%
    filter(type=="unmutated") %>%
    group_by(identifier) %>%
    summarise(unmutated_count=sum(count)) %>%
    ungroup
  mutated <- all_data[[sample_name]] %>%
    filter(type!="unmutated") %>%
    filter(type!="mutation") %>%
    group_by(identifier) %>%
    summarise(mutated_count=sum(count)) %>%
    ungroup
  df <- full_join(mutated,unmutated) %>%
    full_join(shannon_entropy) %>%
    mutate(dev_stage=tmp$dev_stage,sample_name=tmp$`Sample Name`)
  
  data_summary <- rbind(data_summary,
                        df)}
write_tsv(data_summary,"./data/050_data_summary.txt")
