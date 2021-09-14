library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./data/000_SraRunTable.txt") %>%
        select(Run,`Library Name`)

for (run in sra$Run) {
  print(run)
  system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s, 
                -O ./data/010_raw", run))
}
