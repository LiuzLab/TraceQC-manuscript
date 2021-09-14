library(readr)
library(dplyr)
library(stringr)

bulk_samples <- filter(sra,LibrarySource=="GENOMIC")
dir_name <- "./data/010_raw/"

for (i in 1:nrow(bulk_samples)) {
  sample <- bulk_samples$Run[i]
  print(sample)
  system(paste("/mnt/hdd/jasper/Bioinfo_tools/FLASH2-master/flash2",
               sprintf("%s%s_1.fastq",dir_name,sample),
               sprintf("%s%s_2.fastq",dir_name,sample),
               "-o",bulk_samples$Expriment[i],"-d","./data/011_merged_fastq/",
               "-M","250","-O","-m","20","-x","0.1"
  ))
}

for (f in list.files("./data/011_merged_fastq/")) {
  if (str_detect(f,"extendedFrags.fastq")) {
    fname <- strsplit(f,"\\.")[[1]][1]
    system(sprintf("mv ./data/011_merged_fastq/%s ./data/011_merged_fastq/",f,fname))
  }
  else {
    system(sprintf("rm ./data/011_merged_fastq/%s",f))
  }}
