library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./data/000_SraRunTable.txt") %>%
  select(Run,`Library Name`,`Sample Name`,TISSUE,Sampling_replicate)

for (run in sra$Run) {
  print(run)
  system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s 
                -O ./010_raw", run))
}

all_files <- list.files("./data/010_raw/")
file_missing <- c()
for (run in sra$Run) {
  is_downloaded <- sprintf("%s_1.fastq",run) %in% all_files
  if (!is_downloaded) {
    file_missing <- c(file_missing,as.character(run))
    system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s -O ./data/010_raw", run))
  }}
