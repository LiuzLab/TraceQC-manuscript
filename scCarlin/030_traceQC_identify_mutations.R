library(readr)
library(TraceQC)
library(dplyr)
library(tools)

alignment_threshold <- sequence_permutation(ref_file="./data/000_ref.txt")
write_tsv(alignment_threshold,"./data/030.1_alignment_threshold.txt")
alignment_permutation <- read_tsv("./data/030.1_alignment_threshold.txt")
model <- loess(score~permutate_percent,data=alignment_permutation)
alignment_score_cutoff <- predict(model,0.4)

for (fname in list.files("./data/020_alignment_10x/",pattern="Amplicon.txt")) {
  df <- read_tsv(sprintf("./data/020_alignment_10x/%s",fname)) %>%
      filter(CB!="/",UB!="/",score>alignment_score_cutoff)
  
  UMI <- get_read_count_per_UB(df) %>%
      filter(read_count_per_UMI>10)
  df <- inner_join(df,UMI)
  df <- group_by(df,target_seq,target_ref,CB,UB) %>%
    summarise(count=n(),read_count_per_UMI=mean(read_count_per_UMI)) %>%
    ungroup %>%
    filter(count>0.01*read_count_per_UMI)
    
  seq_to_char_input <- group_by(df,target_seq,target_ref) %>%
      summarise() %>%
      ungroup %>%
      mutate(score=1,count=1)
  mutation_event <- seq_to_character(seq_to_char_input,abundance_cutoff=0,
                                   use_CPM=FALSE,alignment_score_cutoff=0)
  df <- select(mutation_event,-alignment_score,-count) %>%
      right_join(df)
  write_tsv(df,sprintf("./data/030.2_mutation_event/%s",fname))}
