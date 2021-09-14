library(readr)
library(dplyr)
library(stringr)

sra <- read_csv("./000_SRA/SraRunTable.txt")
exp <- read_tsv("./000_SRA/Exp_table.txt")

sra <- select(sra,Run,Sample=`Sample Name`,Developmental_stage,LibrarySource,source_name,TISSUE) %>%
  left_join(exp)

for (run in sra$Run) {
  print(run)
  system(sprintf("~/Bioinfotools/sratoolkit.2.10.8-ubuntu64/bin/fasterq-dump %s
                -t ./tmp,  -O ./data/010_raw", run))
}

cell_culture_data <- filter(sra,str_detect(Experiment,"cell_culture")) %>%
  as.data.frame()
dir_name <- "/home/humble_local_25t/jasper/projects/TraceQC-manuscript/GESTALT/data/010_raw/"

# for (i in 1:nrow(cell_culture_data)) {
#   sample <- cell_culture_data$Run[i]
#   print(sample)
#   system(paste("/home/humble_local_25t/jasper/Bioinfotools/FLASH2/flash2",
#                sprintf("%s%s_1.fastq",dir_name,sample),
#                sprintf("%s%s_2.fastq",dir_name,sample),
#                "-o",cell_culture_data$Experiment[i],
#                "-d","/home/humble_local_25t/jasper/projects/TraceQC-manuscript/GESTALT/data/merged_fastq/",
#                "-M","290","-O","-m","20","-x","0.1"
#   ))
# }
