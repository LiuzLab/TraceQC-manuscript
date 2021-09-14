library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

########## aggregate data ##########
sra <- read_csv("./data/000_SraRunTable.txt")
all_data <- list(info_table=sra)

for (run in sra$Run) {
  print(run)
  all_files <- list.files("./data/030.2_traceQC_obj/",pattern=run)
  df <- data.frame()
  for (f in all_files) {
    name <- tools::file_path_sans_ext(f)
    identifier <- strsplit(name,split="_")[[1]][2]
    fname <- list.files("./data/030.3_alignment_threshold/",pattern=sprintf("%s.txt",identifier))
    alignment_permutation <- read_tsv(sprintf("./data/030.3_alignment_threshold/%s",fname))
    model <- loess(score~permutate_percent,data=alignment_permutation)
    
    tmp <- read_tsv(sprintf("./data/030.2_traceQC_obj/%s",f)) %>%
      mutate(identifier=identifier) %>%
      filter(alignment_score>predict(model,0.4))

    df <- rbind(df,tmp)}
  if (nrow(df)>0) {
    all_data[[run]] <- df
  }}
saveRDS(all_data,"./data/040.1_all_data.rds")

########## normalization and filtering ########## 
all_data <- readRDS("./data/040.1_all_data.rds")

all_data_normed <- list(info_table=sra)
tot_read_count_threshold <- 5000
abundance_threhold <- 1e-5

for (run in names(all_data)[2:length(all_data)]) {
  df <- all_data[[run]]
  df_tmp <- group_by(df,target_seq,identifier) %>%
    summarise(count=mean(count))
  if (sum(df_tmp$count)<5000) next
  print(paste(run,sum(df_tmp$count)))
  df_tmp$count <- df_tmp$count * 1e6 / sum(df_tmp$count)
  threshold <- 1e6 * abundance_threhold
  df_tmp <- filter(df_tmp,count>threshold)
  
  df <- select(df,-count) %>%
    right_join(df_tmp)
  all_data_normed[[run]] <- df
}
saveRDS(all_data_normed,"./data/040.2_all_data_normed.rds")

