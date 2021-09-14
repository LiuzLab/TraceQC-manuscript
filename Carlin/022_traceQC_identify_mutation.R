library(TraceQC)
library(ggplot2)

alignment_permutation <- read_tsv("./data/021_alignment_threshold.txt")
model <- loess(score~permutate_percent,data=alignment_permutation)
alignment_score_cutoff <- predict(model,0.4)
for (input_file in list.files("./data/020_alignment/")) {
  name <- tools::file_path_sans_ext(input_file)
  aligned_reads <- read_tsv(sprintf("./data/020_alignment/%s",input_file)) %>%
    filter(score>alignment_score_cutoff) %>%
    group_by(target_seq,target_ref) %>%
    summarise(count=n(),score=max(score)) %>%
    ungroup
  out_file <- paste("./data/022_traceqQC_obj/",name,".txt",sep="")
  print(name)
  mutation_event <- seq_to_character(aligned_reads,abundance_cutoff=1e-5,
                        use_CPM=TRUE,alignment_score_cutoff=alignment_score_cutoff)
  write_tsv(mutation_event,out_file)}
