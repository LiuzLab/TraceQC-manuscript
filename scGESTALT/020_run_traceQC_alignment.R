library(TraceQC)
library(fastqcr)
library(tools)
library(stringr)
library(readr)
library(ggpubr)
library(dplyr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt") %>%
  select(Run,Developmental_stage,Experiment,`GEO_Accession (exp)`,
         LibraryLayout,LibrarySource,LibrarySource,source_name)
exp <- read_tsv("./data/000_SRA/Exp_table.txt") %>%
  select(`GEO_Accession (exp)`=Sample,Exp=Experiment)
sra <- left_join(sra,exp)

ref_file <- "./data/000_ref.txt"
ref <- parse_ref_file(ref_file)
plot_construct(ref)
runs <- c("SRR6176748","SRR6176749","SRR6176750")

for (run in runs) {
    out_alignment_file <- paste("./data/020.1_alignment/",run,".txt",sep="")
    sequence_alignment(input_file=paste("./data/001_raw_divided/single_end/",run,".fastq",sep=""),
              ref_file=ref_file,output_file=out_alignment_file,ncores=16)}
