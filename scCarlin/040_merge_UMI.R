library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(TraceQC)

for (fname in list.files("./data/030.2_mutation_event/")) {
  print(fname)
  df <- read_tsv(sprintf("./data/030.2_mutation_event/%s",fname)) %>%
    group_by(CB,UB,type,start,length,mutate_to) %>%
    summarise(count=sum(count),read_count_per_UMI=mean(read_count_per_UMI)) %>%
    ungroup
  
  df <- group_by(df,CB,UB) %>%
    nest() %>%
    mutate(mutations=map(data,filter_mutations,TRUE,0.5)) %>%
    unnest(mutations) %>%
    select(-data)
  
  UMI <- get_UMI_count_per_CB(df)
  
  df <- group_by(df,CB,type,start,length,mutate_to) %>%
    summarise(UMI_count=n()) %>%
    ungroup %>%
    left_join(UMI) %>%
    filter(UMI_count >= 0.5*UMI_per_CB)
  write_tsv(df,sprintf("./data/040_cell_mutations/%s",fname))
}
