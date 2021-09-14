library(readr)
library(TraceQC)
library(dplyr)
library(tools)

alignment_threshold <- sequence_permutation(ref_file="./data/000_ref.txt",
                                      penalize_end_gaps=0,read_length=100)
write_tsv(alignment_threshold,"./data/031_alignment_threshold.txt")
model <- loess(score~permutate_percent,data=alignment_threshold)
alignment_score_cutoff <- predict(model,0.4)

for (fname in list.files("./data/030_traceQC_alignment/")) {
    df <- read_tsv(sprintf("./data/030_traceQC_alignment/%s",fname))
    if (nrow(df)>0) {
      df <- filter(df,CB!="/",UB!="/",score>alignment_score_cutoff)
  
    UMI <- get_read_count_per_UB(df)
    df <- inner_join(df,UMI)
    df <- group_by(df,target_seq,target_ref,CB,UB) %>%
      summarise(count=n(),read_count_per_UMI=mean(read_count_per_UMI)) %>%
      ungroup
    
    seq_to_char_input <- group_by(df,target_seq,target_ref) %>%
      summarise() %>%
      ungroup %>%
      mutate(score=1,count=1)
    mutation_event <- seq_to_character(seq_to_char_input,abundance_cutoff=0,
                                       use_CPM=FALSE,alignment_score_cutoff=0)
    df <- select(mutation_event,-alignment_score,-count) %>%
      right_join(df)
    write_tsv(df,sprintf("./data/040_mutation_event/%s",fname))}}
