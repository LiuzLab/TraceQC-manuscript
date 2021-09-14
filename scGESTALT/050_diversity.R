library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)
library(tools)

sra <- read_csv("./data/000_SRA/SraRunTable.txt") %>%
  select(Run,Developmental_stage,Experiment,`GEO_Accession (exp)`,
         LibraryLayout,LibrarySource,LibrarySource,source_name)
exp <- read_tsv("./data/000_SRA/Exp_table.txt") %>%
  select(`GEO_Accession (exp)`=Sample,Exp=Experiment)
sra <- left_join(sra,exp)

all_data <- list(info_table=sra)
for (input_file in list.files("./data/040_cell_mutations/")) {
  run <- file_path_sans_ext(input_file)
  print(run)
  all_data[[run]] <- read_tsv(sprintf("./data/040_cell_mutations/%s",input_file))
}
saveRDS(all_data,"./data/050.1_all_data.rds")
all_data <- readRDS("./data/050.1_all_data.rds")

data_summary <- data.frame()
cell_barcodes <- read_delim("./data/002_cell_barcodes.txt",delim=" ") %>%
  separate(InDrops_cellBarcode,c("tissue","CB"),sep="_") %>%
  select(-Shortname)
unique_CB <- (table(cell_barcodes$CB) == 1)
cell_barcodes <- filter(cell_barcodes,unique_CB[CB])

for (run in names(all_data)[2:length(all_data)]) {
  myrun <- all_data[[run]]
  for (mytissue in unique(myrun$tissue)) {
    sample <- filter(myrun,tissue==mytissue)
    tmp <- sample %>%
      mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
      group_by(cigar) %>%
      summarise(count=n())
    counts <- tmp$count
    
    number_of_cells <- length(unique(sample$CB))
    umi <- group_by(sample,CB) %>%
      summarise(umi_count=mean(tot_num_of_UMI)) %>%
      ungroup
    shannon_entropy <- Entropy(counts)
    gini_index <- Gini(counts)
    diversity <- length(counts)
    
    transcriptome_barcodes <- filter(cell_barcodes,tissue==mytissue)
    
    data_summary <- rbind(data_summary,
                          data.frame(run=run,tissue=mytissue,
                                     number_of_cells=number_of_cells,
                                     umi_per_cell=mean(umi$umi_count),
                                     shannon_entropy=shannon_entropy,
                                     gini_index=gini_index,
                                     diversity=diversity,
                                     number_of_cells_transcriptome=length(transcriptome_barcodes$CB),
                                     cells_overlap=length(intersect(transcriptome_barcodes$CB,umi$CB))
                          ))}}
write_tsv(data_summary,"./data/050.2_data_summary.txt")
data_summary <- read_tsv("./data/050.2_data_summary.txt")
