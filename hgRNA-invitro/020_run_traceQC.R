library(readr)
library(dplyr)
library(stringr)
library(tools)
library(TraceQC)

sra <- read_csv("./data/000_SraRunTable.txt") %>%
  select(Run,`Library Name`)
fastq_dir <- "./data/010_raw/"

all_refs <- file_path_sans_ext(list.files("./data/000_ref/"))

for (i in 1:nrow(sra)) {
  prefix <- strsplit(as.character(sra[i,"Library Name"]),split="-")[[1]][1]
  if (prefix %in% all_refs) {
    input_file <- sprintf("./data/010_raw/%s.fastq",sra[i,"Run"])
    ref_file <- sprintf("./data/000_ref/%s.txt",prefix)
    output_file <- sprintf("./data/020_alignment/%s.txt",sra[i,"Library Name"])
    print(paste(ref_file,output_file))
    sequence_alignment(input_file=input_file,ref_file=ref_file,
                       output_file=output_file)
  }}
