library(readr)
library(dplyr)
library(stringr)
library(TraceQC)
library(fastqcr)
library(tools)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Sample=`GEO_Accession (exp)`,Developmental_stage,LibrarySource,source_name,TISSUE) %>%
  left_join(exp)

################ cell_culture ################ 
cell_culture_data <- filter(sra,str_detect(Experiment,"cell_culture")) %>%
  as.data.frame()

ref_file <- "./data/000_ref_v1.txt"
qc_dir <- "./fastqc/"
fastq_dir <- "./data/010_raw/"
fastqc("./data/010_raw/",qc.dir=qc_dir)

for (run in c("SRR3561135",cell_culture_data$Run)) {
  out_alignment_file <- paste("./data/020.1_alignment/",run,".txt",sep="")
  if (!file.exists(out_alignment_file)) {
    print(run)
    sequence_alignment(input_file=sprintf("%s%s_1.fastq",fastq_dir,run),
                       ref_file=ref_file,output_file=out_alignment_file,ncores=16)}
}

alignment_threshold <- sequnce_alignment_threshold(ref_file=ref_file)
write_tsv(alignment_threshold,"./data/020.2_alignment_threshold_cell_culture.txt")
alignment_score_cutoff <- filter(alignment_threshold,permutate_percent==0.4) %>%
  pull(score) %>% min()

for (run in c("SRR3561135",cell_culture_data$Run)) {
  aligned_reads <- read_tsv(sprintf("./data/020.1_alignment/%s.txt",run)) %>%
    filter(score>alignment_score_cutoff) %>%
    group_by(target_seq,target_ref) %>%
    summarise(count=n(),score=max(score)) %>%
    ungroup
  mutation_event <- seq_to_character(aligned_reads,ncores=1,
                                     use_CPM=TRUE,alignment_score_cutoff=alignment_score_cutoff,
                                     abundance_cutoff=1e-5)
  write_tsv(mutation_event,sprintf("./data/020.3_mutation_event/cell_culture/%s.txt",run))}

################ in-vivo ################ 
whole_organism_data <- filter(sra,source_name=="whole_organism"&
                                !str_detect(Experiment,"v6"))

fastq_dir <- "./data/010_raw/"
ref <- parse_ref_file("./data/000_ref_v6.txt")
plot_construct(ref)
ref_r1 <- "./data/000_ref_v6_r1.txt"
ref_r2 <- "./data/000_ref_v6_r2.txt"
for (run in whole_organism_data$Run) {
  print(run)
  out_alignment_file_r1 <- paste("./data/020.1_alignment/",run,"_1.txt",sep="")
  out_alignment_file_r2 <- paste("./data/020.1_alignment/",run,"_2.txt",sep="")
  sequence_alignment(input_file=sprintf("%s%s_1.fastq",fastq_dir,run),
                     ref_file=ref_r1,output_file=out_alignment_file_r1,ncores=10)
  sequence_alignment(input_file=sprintf("%s%s_2.fastq",fastq_dir,run),
                     ref_file=ref_r2,output_file=out_alignment_file_r2,ncores=10)
}

alignment_threshold <- sequence_permutation(ref_file=ref_r1)
model <- loess(score~permutate_percent,data=alignment_threshold)
alignment_score_cutoff_r1 <- predict(model,0.4)
write_tsv(alignment_threshold,"./data/020.2_alignment_threshold_r1.txt")
alignment_threshold <- sequence_permutation(ref_file=ref_r1)
model <- loess(score~permutate_percent,data=alignment_threshold)
alignment_score_cutoff_r2 <- predict(model,0.4)
write_tsv(alignment_threshold,"./data/020.2_alignment_threshold_r2.txt")

for (run in whole_organism_data$Run) {
  print(run)
  aligned_reads_r1 <- read_tsv(sprintf("./data/020.1_alignment/%s_1.txt",run)) %>%
    filter(score>alignment_score_cutoff_r1) %>%
    group_by(target_seq,target_ref) %>%
    summarise(count=n(),score=max(score)) %>%
    ungroup
  mutation_event <- seq_to_character(aligned_reads_r1,
                                     use_CPM=TRUE,alignment_score_cutoff=alignment_score_cutoff_r1,
                                     abundance_cutoff=1e-5)
  write_tsv(mutation_event,sprintf("./data/020.3_mutation_event/whole_organism/%s_1.txt",run))
  
  aligned_reads_r2 <- read_tsv(sprintf("./data/020.1_alignment/%s_2.txt",run)) %>%
    filter(score>alignment_score_cutoff_r2) %>%
    group_by(target_seq,target_ref) %>%
    summarise(count=n(),score=max(score)) %>%
    ungroup
  mutation_event <- seq_to_character(aligned_reads_r2,
                                     use_CPM=FALSE,alignment_score_cutoff=alignment_score_cutoff_r2,
                                     abundance_cutoff=1e-5)
  write_tsv(mutation_event,sprintf("./data/020.3_mutation_event/whole_organism/%s_2.txt",run))}
