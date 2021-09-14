# library(readr)
# library(tidyr)
# 
cell_barcodes <- read_delim("./data/002_cell_barcodes.txt",delim=" ") %>%
  separate(InDrops_cellBarcode,c("tissue","CB"),sep="_") %>%
  select(-Shortname)
unique_CB <- (table(cell_barcodes$CB) == 1)
cell_barcodes <- filter(cell_barcodes,unique_CB[CB])
# 
# for (fname in list.files("./data/030_mutation_event/")) {
#   print(fname)
#   df <- read_tsv(sprintf("./data/030_mutation_event/%s",fname)) %>%
#     separate(CB,c("CB","suffix"),sep=":") %>%
#     inner_join(cell_barcodes) %>%
#     select(-suffix) %>%
#     group_by(CB,UB,tissue,type,start,length,mutate_to) %>%
#     summarise() %>%
#     ungroup
#   
#   UMI <- group_by(df,CB,UB) %>%
#     summarise() %>%
#     ungroup %>%
#     group_by(CB) %>%
#     summarise(tot_num_of_UMI=n()) %>%
#     ungroup
#   
#   df <- group_by(df,CB,tissue,type,start,length,mutate_to) %>%
#     summarise(UMI_count=n()) %>%
#     ungroup %>%
#     left_join(UMI) %>%
#     filter(UMI_count >= 0.5*tot_num_of_UMI)
#   write_tsv(df,sprintf("./data/040_cell_mutations/%s",fname))
# }

library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(TraceQC)

for (fname in list.files("./data/030_mutation_event/")) {
  print(fname)
  df <- read_tsv(sprintf("./data/030_mutation_event/%s",fname)) %>%
    separate(CB,c("CB","suffix"),sep=":") %>%
    inner_join(cell_barcodes) %>%
    select(-suffix) %>%
    group_by(CB,UB,tissue,type,start,length,mutate_to) %>%
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




