library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)
library(tools)

sra <- read_csv("./data/000_SRA/SraRunTable.txt") %>%
  select(Run,Developmental_stage,Sample=`GEO_Accession (exp)`,LibrarySource,source_name,TISSUE,Time_point)

exp <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- left_join(sra,exp)

################ in-vitro ################ 
all_data <- list(info_table=sra)
for (input_file in list.files("./data/050_cell_mutations/")) {
  run <- file_path_sans_ext(input_file)
  print(run)
  all_data[[run]] <- read_tsv(sprintf("./data/050_cell_mutations/%s",input_file))
}
saveRDS(all_data,"./data/060.1_all_data.rds")

######################################
data_summary <- data.frame()
for (run in names(all_data)[2:length(all_data)]) {
  experiment <- filter(sra,Run==run) %>% pull(Experiment)
  sample_name <- experiment
  
  sample <- all_data[[run]] %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(cigar) %>%
    summarise(count=n())
  counts <- sample$count
  
  number_of_cells <- length(unique(all_data[[run]]$CB))
  shannon_entropy <- Entropy(counts)
  gini_index <- Gini(counts)
  diversity <- length(counts)
  transcriptome_barcodes <- read_tsv(sprintf("./data/020_cellranger_output/%s/outs/filtered_feature_bc_matrix/barcodes.tsv",
                                             run),col_names=FALSE)
  
  data_summary <- rbind(data_summary,
                        data.frame(sample_name=sample_name,
                                   shannon_entropy=shannon_entropy,
                                   gini_index=gini_index,
                                   diversity=diversity,
                                   number_of_cells=number_of_cells,
                                   number_of_cells_transcriptome=length(transcriptome_barcodes$X1),
                                   cells_overlap=length(intersect(transcriptome_barcodes$X1,
                                                                  unique(all_data[[run]]$CB)))))}
write_tsv(data_summary,"./data/060.2_data_summary.txt")
data_summary <- read_tsv("./data/060.2_data_summary.txt")
