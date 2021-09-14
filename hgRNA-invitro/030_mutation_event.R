library(TraceQC)
library(tools)

sra <- read_csv("./data/000_SraRunTable.txt") %>%
  select(Run,`Library Name`)

for (ref_file in list.files("./data/000_ref/")) {
  alignment_threshold <- sequence_permutation(ref_file=sprintf("./data/000_ref/%s",ref_file))
  write_tsv(alignment_threshold,sprintf("./data/030.1_alignment_threshold/alignment_threshold_%s",ref_file))}

for (input_file in list.files("./data/020_alignment/")) {
  name <- file_path_sans_ext(input_file)
  prefix <- strsplit(name,split="-")[[1]][1]
  ref_file <- sprintf("./data/000_ref/%s.txt",prefix)
  alignment_threshold <- read_tsv(sprintf("./data/030.1_alignment_threshold/alignment_threshold_%s.txt",prefix)) %>%
    filter(permutate_percent==0.4) %>%
    pull(score) %>%
    min()
    
  fastqc_file <- sprintf("./fastqc/%s_fastqc.zip",
                         filter(sra,`Library Name`==name) %>% pull(Run))
  out_file <- sprintf("./data/030.2_traceQC_obj/%s.rds",name)
  obj <- create_TraceQC_object(aligned_reads_file = paste("./data/020_alignment/",input_file,sep=""),
                               ref_file = ref_file,
                               fastqc_file = fastqc_file,
                               use_CPM=TRUE,
                               alignment_score_cutoff=alignment_threshold,
                               abundance_cutoff=1e-5,
                               ncores = 1)
  saveRDS(obj,file=out_file)}
