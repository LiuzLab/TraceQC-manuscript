library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")

for (run in sra$Run) {
  print(run)
  system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s -O ./data/010_raw", run))
}

# files <- list.files("./raw/")
# files_missing <- c()
# for (i in 1:nrow(sra)) {
#   if (sra[i,"LibraryLayout"]=="SINGLE") {
#     is_downloaded <- sprintf("%s.fastq",sra[i,"Run"]) %in% files
#   } else if (sra[i,"LibraryLayout"]=="PAIRED") {
#     is_downloaded <- sprintf("%s_1.fastq",sra[i,"Run"]) %in% files
#   }
#   if (!is_downloaded) {
#     files_missing <-(sra[i,"Run"]))
#   # system(sprintf("~/hdd/Bioinfo_tools/sratoolkit.2.10.0-ubuntu64/bin/prefetch %s && ~/hdd/Bioinfo_tools/sratoolkit.2.10.0-ubuntu64/bin/fasterq-dump %s",sra[i,"Run"],sra[i,"Run"]))
#   }
# }
