library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)
library(tools)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp_table <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Cell_type,Sample=`Sample Name`,LibrarySource,TREATMENT) %>%
  left_join(exp_table)

################ in-vitro ################ 
all_data <- list(info_table=sra)
for (input_file in list.files("./data/040_cell_mutations/")) {
  run <- file_path_sans_ext(input_file)
  print(run)
  all_data[[run]] <- read_tsv(sprintf("./data/040_cell_mutations/%s",input_file))
}
saveRDS(all_data,"./data/050.1_all_data.rds")
all_data <- readRDS("./data/050.1_all_data.rds")

######################################
for (sample_name in list.files("./data/010_cellranger_output/")) {
  tsv_file <- sprintf("./data/010_cellranger_output/%s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
                      sample_name)
  system(sprintf("gzip -d %s",tsv_file))}
