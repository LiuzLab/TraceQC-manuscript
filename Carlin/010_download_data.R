library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp <- read_tsv("./data/000_SRA/Exp_table.txt")

sra <- select(sra,Run,Cell_type,Sample=`Sample Name`,LibrarySource,TREATMENT) %>%
  left_join(exp)

for (run in sra$Run) {
  print(run)
  system(sprintf("~/hdd/Bioinfo_tools/sratoolkit.2.10.0-ubuntu64/bin/fasterq-dump %s -O ./data/010_raw", run))
}

all_files <- list.files("./data/010_raw/")
file_missing <- c()
for (run in sra$Run) {
  is_downloaded <- sprintf("%s_1.fastq",run) %in% all_files
  if (!is_downloaded) {
    file_missing <- c(file_missing,as.character(run))
    system(sprintf("~/hdd/Bioinfo_tools/sratoolkit.2.10.0-ubuntu64/bin/fasterq-dump %s -O ./data/010_raw", run))
  }}
