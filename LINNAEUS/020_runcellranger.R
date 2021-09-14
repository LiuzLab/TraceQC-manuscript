library(readr)
library(dplyr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt") %>%
  select(Run,Developmental_stage,`GEO_Accession (exp)`,LibrarySource,source_name,TISSUE,Time_point)

exp_table <- read_tsv("./data/000_SRA/Exp_table.txt")

prefix <- "/home/humble_local_25t/jasper/projects/TraceQC-manuscript/LINNAEUS"
setwd(paste(prefix,"/data/020_cellranger_output",sep=""))
sc_10x <- filter(sra,LibrarySource=="TRANSCRIPTOMIC")

for (run in unique(sc_10x$Run)) {
  system(sprintf("mkdir -p %s/data/011_cellranger_input/%s",prefix,run))
  system(sprintf("ln -rs %s/data/010_raw/%s_1.fastq %s/data/011_cellranger_input/%s/%s_S1_L001_R1_001.fastq",prefix,run,prefix,run,run))
  system(sprintf("ln -rs %s/data/010_raw/%s_2.fastq %s/data/011_cellranger_input/%s/%s_S1_L001_R2_001.fastq",prefix,run,prefix,run,run))
  
  cmd <- paste("/home/humble_local_25t/jasper/Bioinfotools/cellranger-3.1.0/cellranger count",
               sprintf("--id=%s --fastqs=%s/data/011_cellranger_input/%s/",run,prefix,run),
               "--transcriptome=/home/humble_local_25t/jasper/Bioinfotools/cellranger-3.1.0/ref/Danio.rerio_genome",
               sprintf("--sample=%s",run))
  system(cmd)
}
