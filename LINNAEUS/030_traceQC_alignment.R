library(TraceQC)
library(fastqcr)
library(tools)
library(dplyr)

ref_file <- "./data/000_ref.txt"
for (run in list.files("./data/021_cellranger_bam_filtered/")) {
  print(run)
  input_file <- sprintf("./data/021_cellranger_bam_filtered/%s",run)
  out_alignment_file <- paste("./data/030_traceQC_alignment/",file_path_sans_ext(run),".txt",sep="")
  sequence_alignment_for_10x(input_file,ref_file,penalize_end_gaps=0,output_file=out_alignment_file,ncores=16)}
