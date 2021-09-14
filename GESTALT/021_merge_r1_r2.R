library(readr)
library(dplyr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Sample=`Sample Name`,Developmental_stage,LibrarySource,source_name,TISSUE) %>%
  left_join(exp)
whole_organism_data <- filter(sra,source_name=="whole_organism"&
                                !str_detect(Experiment,"v6"))

for (run in whole_organism_data$Run) {
  r1 <- read_tsv(sprintf("./data/020.3_mutation_event/whole_organism/%s_1.txt",run)) %>%
    mutate(read="r1")
  
  r2 <- read_tsv(sprintf("./data/020.3_mutation_event/whole_organism/%s_2.txt",run)) %>%
    mutate(start=274-start-length,read="r2")
  
  df <- rbind(r1,r2)
  write_tsv(df,sprintf("./data/020.3_mutation_event/whole_organism_merged/%s.txt",run))}
