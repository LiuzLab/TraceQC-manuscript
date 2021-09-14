library(readr)
library(dplyr)
library(TraceQC)

# for (sample_name in list.files("./data/020_cellranger_output/")) {
#   tsv_file <- sprintf("./data/020_cellranger_output/%s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#                       sample_name)
#   system(sprintf("gzip -d %s",tsv_file))}

ref <- parse_ref_file("./data/000_ref.txt")
target_start <- filter(ref$regions,region=="target") %>% pull(start)
target_end <- filter(ref$regions,region=="target") %>% pull(end)
sgRNA_start <- filter(ref$regions,region=="sgRNA") %>% pull(start)
sgRNA_end <- filter(ref$regions,region=="sgRNA") %>% pull(end)

for (fname in list.files("./data/040_mutation_event/")) {
  print(fname)
  df <- read_tsv(sprintf("./data/040_mutation_event/%s",fname)) %>%
    filter(!(type=="deletion"&(start==1|(start+length)==target_end-target_start+2))) %>%
    filter(start>(sgRNA_start-20)) %>%
    group_by(CB,UB,type,start,length,mutate_to) %>%
    summarise(count=sum(count),read_count_per_UMI=mean(read_count_per_UMI)) %>%
    ungroup
  # filtered_cell_barcodes <- read_tsv(sprintf("./data/020_cellranger_output/%s/outs/filtered_feature_bc_matrix/barcodes.tsv",
  #                                   tools::file_path_sans_ext(fname)),col_names=FALSE) %>%
  #   pull(X1)
  # df <- filter(df,CB %in% filtered_cell_barcodes)
  UMI <- get_UMI_count_per_CB(df)
  
  df <- group_by(df,CB,type,start,length,mutate_to) %>%
    summarise(UMI_count=n()) %>%
    ungroup %>%
    left_join(UMI)
  write_tsv(df,sprintf("./data/050_cell_mutations/%s",fname))}
