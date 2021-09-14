library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

############## calculate measures ############## 
all_data <- readRDS("./data/030_all_data.rds")
data_summary <- data.frame()
for (sample_name in names(all_data)[2:length(all_data)]) {
    sample <- all_data[[sample_name]] %>%
      mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
      group_by(cigar) %>%
      summarise(count=sum(count))
    counts <- sample$count
    
    tmp <- filter(all_data[[sample_name]],type!="unmutated") %>%
      group_by(target_seq,type) %>%
      summarise(length=sum(length),count=mean(count)) %>%
      ungroup %>%
      spread(type,length) %>%
      mutate(deletion=replace_na(deletion,0),
             insertion=replace_na(insertion,0),
             substitution=replace_na(substitution,0))
    avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
    
    shannon_entropy <- Entropy(counts)
    gini_index <- Gini(counts)
    diversity <- length(counts)
    unmutated_count <- all_data[[sample_name]] %>%
      filter(type=="unmutated") %>%
      pull(count)
    data_summary <- rbind(data_summary,
                data.frame(sample_name=sample_name,
                           avg_deletion_length=avg_deletion_length,
                           shannon_entropy=shannon_entropy,
                           gini_index=gini_index,
                           diversity=diversity,
                           unmutated_percentage=unmutated_count/1e6))}
write_tsv(data_summary,"./data/040_data_summary.txt")
