library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt") %>%
  select(Run,Developmental_stage,Experiment,`GEO_Accession (exp)`,
         LibraryLayout,LibrarySource,LibrarySource,source_name)
exp <- read_tsv("./data/000_SRA/Exp_table.txt") %>%
  select(`GEO_Accession (exp)`=Sample,Exp=Experiment)
sra <- left_join(sra,exp)

paired <- filter(sra,LibraryLayout=="PAIRED")
for (run in paired$Run) {
  print(run)
  system(sprintf("~/hdd/Bioinfo_tools/sratoolkit.2.10.0-ubuntu64/bin/fasterq-dump %s -O ./data/000_raw/", run))
}

single <- filter(sra,LibraryLayout=="SINGLE")
all_files <- list.files("./data/000_raw/")
for (run in single$Run) {
  if (!(sprintf("%s.fastq",run) %in% all_files)){
    print(run)
    # system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/prefetch %s && ~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s", run, run))
  }
}

for (run in paired$Run) {
  system(sprintf("mv ./data/000_raw/%s_1.fastq ./data/001_raw_divided/paired_end/%s_R1.fastq", run, run))
  system(sprintf("mv ./data/000_raw/%s_2.fastq ./data/001_raw_divided/paired_end/%s_R2.fastq", run, run))
}

for (run in single$Run) {
  system(sprintf("mv ./data/000_raw/%s.fastq ./data/001_raw_divided/single_end/%s.fastq", run, run))
}
