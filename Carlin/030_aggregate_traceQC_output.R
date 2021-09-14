library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

############## read traceQC output ############## 
sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Cell_type,Sample=`Sample Name`,LibrarySource,TREATMENT) %>%
  left_join(exp)

all_data <- list(info_table=sra)
for (f in list.files("./data/022_traceqQC_obj/")) {
  if (!(str_detect(f,"pop"))) {
  name <- tools::file_path_sans_ext(f)
  print(name)
  all_data[[name]] <- read_tsv(sprintf("./data/022_traceqQC_obj/%s",f))}
}
saveRDS(all_data,"./data/030_all_data.rds")
