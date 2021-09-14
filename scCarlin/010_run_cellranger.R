library(readr)
library(dplyr)

system("ln -rs ../Carlin/data/010_raw/* ./data/000_raw/")

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp_table <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Cell_type,Sample=`Sample Name`,LibrarySource,TREATMENT) %>%
  left_join(exp_table)

prefix <- "/home/humble_local_25t/jasper/projects/TraceQC-manuscript/scCarlin/"
sc_10x <- filter(sra,LibrarySource=="TRANSCRIPTOMIC",Cell_type=="Bone marrow LSK")

for (experiment in unique(sc_10x$Experiment)) {
  samples <- filter(sc_10x,Experiment==experiment)
  system(sprintf("mkdir -p ./data/001_cellranger_input/%s",experiment))
  for (i in 1:nrow(samples)) {
    run <- samples[i,"Run"]
    system(sprintf("ln -rs ./data/000_raw/%s_1.fastq ./data/001_cellranger_input/%s/%s_S1_L00%d_R1_001.fastq",run,experiment,experiment,i))
    system(sprintf("ln -rs ./data/000_raw/%s_2.fastq ./data/001_cellranger_input/%s/%s_S1_L00%d_R2_001.fastq",run,experiment,experiment,i))
  }
  
  cmd <- paste("/home/humble_local_25t/jasper/Bioinfotools/cellranger-3.1.0/cellrangers count",
               sprintf("--id=%s --fastqs=~/hdd/TraceQC-manuscript/scCarlin/data/001_cellranger_input/%s/",experiment,experiment),
               "--transcriptome=/mnt/hdd/jasper/Bioinfo_tools/cell_ranger/ref/refdata-cellranger-mm10-3.0.0",
               sprintf("--sample=%s",experiment))
  system(cmd)}
