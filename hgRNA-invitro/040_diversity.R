library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

############## read traceQC output ############## 
sra <- read_csv("./data/000_SraRunTable.txt")
all_data <- list(info_table=sra)
for (f in list.files("./data/030.2_traceQC_obj/")) {
  name <- tools::file_path_sans_ext(f)
  print(name)
  all_data[[name]] <- readRDS(sprintf("./data/030.2_traceQC_obj/%s",f))[["mutation"]]}
saveRDS(all_data,"./data/040.1_all_data.rds")
all_data <- readRDS("./data/040.1_all_data.rds")

############## calculate measures ############## 
data_summary <- data.frame()
for (sample_name in names(all_data)[2:length(all_data)]) {
  tmp <- strsplit(sample_name,"-") %>% unlist()
  barcode <- tmp[1]
  day <- tmp[2]
  if (str_detect(sample_name,"\\(A'\\)")) {
    barcode <- "A'21"
    day <- substr(day,start=1,stop=nchar(day)-4)
  }

  sample <- all_data[[sample_name]] %>%
    filter(type!="mutation") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(cigar) %>%
    summarise(count=sum(count))
  counts <- sample$count
  
  tmp <- filter(all_data[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  
  shannon_entropy <- Entropy(counts)
  gini_index <- Gini(counts)
  diversity <- length(counts)
  unmutated_count <- all_data[[sample_name]] %>%
    filter(type=="unmutated") %>%
    pull(count)
  if (length(unmutated_count)==0) {unmutated_count<-0}
  
  data_summary <- rbind(data_summary,
              data.frame(barcode=barcode,time=day,
                         shannon_entropy=shannon_entropy,
                         avg_deletion_length=avg_deletion_length,
                         gini_index=gini_index,
                         diversity=diversity,
                         unmutated_percentage=unmutated_count/1e6))}
write_tsv(data_summary,"./data/040.2_data_summary.txt")

