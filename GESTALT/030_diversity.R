library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DescTools)

sra <- read_csv("./data/000_SRA/SraRunTable.txt")
exp <- read_tsv("./data/000_SRA/Exp_table.txt")
sra <- select(sra,Run,Sample=`Sample Name`,Developmental_stage,LibrarySource,source_name,TISSUE) %>%
  left_join(exp)

################ in-vitro ################ 
cell_culture_data <- filter(sra,str_detect(Experiment,"cell_culture")) %>%
  as.data.frame()

all_data <- list(info_table=sra)
for (run in cell_culture_data$Run) {
  print(run)
  all_data[[run]] <- read_tsv(sprintf("./data/020.3_mutation_event/cell_culture/%s.txt",run))
}
saveRDS(all_data,"./data/030.1_all_data_cell_culture.rds")
all_data <- readRDS("./data/030.1_all_data_cell_culture.rds")

######################################
data_summary <- data.frame()
for (run in names(all_data)[2:length(all_data)]) {
    experiment <- filter(cell_culture_data,Run==run) %>% pull(Experiment)
    tmp <- strsplit(experiment,split="_") %>% unlist()
    sample_name <- tmp[length(tmp)]
    
    sample <- all_data[[run]] %>%
      mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
      group_by(cigar) %>%
      summarise(count=sum(count))
    counts <- sample$count

    shannon_entropy <- Entropy(counts)
    gini_index <- Gini(counts)
    diversity <- length(counts)
    unmutated_count <- all_data[[run]] %>%
      filter(type=="unmutated") %>%
      pull(count)
    data_summary <- rbind(data_summary,
                data.frame(sample_name=sample_name,
                           shannon_entropy=shannon_entropy,
                           gini_index=gini_index,
                           diversity=diversity,
                           unmutated_percentage=unmutated_count/1e6))}
write_tsv(data_summary,"./data/030.2_invitro_data_summary.txt")
data_summary <- read_tsv("./data/030.2_invitro_data_summary.txt")

######################################
df_tmp <- filter(data_summary,!str_detect(sample_name,"a")&
                                !str_detect(sample_name,"b"))
data_summary <- filter(data_summary,
                str_detect(sample_name,"a")|str_detect(sample_name,"b")) %>%
  mutate(group=sample_name,
         time = 2)

data_summary <- rbind(data_summary,
              mutate(df_tmp,group=paste(sample_name,"a",sep=""),time=1)) %>%
  rbind(mutate(df_tmp,group=paste(sample_name,"b",sep=""),time=1))


ggplot(data=data_summary,aes(x=time,y=unmutated_percentage,
                   group=group,color=group)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=c(1,2)) +
  theme_classic()
ggsave("./figures/invitro_unmutated.png",height=5,width=8)

ggplot(data=data_summary,aes(x=time,y=shannon_entropy,
                   group=group,color=group)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=c(1,2)) +
  theme_classic()
ggsave("./figures/invitro_shannon_entropy.png",height=5,width=8)

ggplot(data=data_summary,aes(x=time,y=gini_index,
                   group=group,color=group)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=c(1,2)) +
  theme_classic()
ggsave("./figures/invitro_gini_index.png",height=5,width=8)

ggplot(data=data_summary,aes(x=time,y=diversity,
                   group=group,color=group)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=c(1,2)) +
  theme_classic()
ggsave("./figures/invitro_diversity.png",height=5,width=8)

################ in-vivo ################
whole_organism_data <- filter(sra,source_name=="whole_organism"&
                                !str_detect(Experiment,"v6"))

all_data <- list(info_table=sra)
for (run in whole_organism_data$Run) {
  print(run)
  all_data[[run]] <- read_tsv(sprintf("./data/020.3_mutation_event/whole_organism_merged/%s.txt",run))
}
saveRDS(all_data,"./data/030.1_all_data_whole_organism.rds")
all_data <- readRDS("./data/030.1_all_data_whole_organism.rds")

######################################
data_summary <- data.frame()
for (i in 1:nrow(whole_organism_data)) {
  cas9_RNP <- strsplit(as.character(whole_organism_data[i,"Experiment"]),"_")[[1]][3]
  dev_stage <- as.character(whole_organism_data[i,"Developmental_stage"])
  run <- as.character(whole_organism_data[i,"Run"])
  
  df <- read_tsv(sprintf("./data/020.3_mutation_event/whole_organism_merged/%s.txt",run))
  sample <- df %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(cigar) %>%
    summarise(count=sum(count))
  counts <- sample$count
  
  shannon_entropy <- Entropy(counts)
  gini_index <- Gini(counts)
  diversity <- length(counts)
  unmutated_count <- df %>%
    filter(type=="unmutated") %>%
    pull(count) %>%
    sum()
  data_summary <- rbind(data_summary,
                        data.frame(run=run,cas9_RNP=cas9_RNP,
                                   dev_stage=dev_stage,
                                   shannon_entropy=shannon_entropy,
                                   gini_index=gini_index,
                                   diversity=diversity,
                                   unmutated_percentage=unmutated_count/1e6))}
write_tsv(data_summary,"./data/030.2_invivo_data_summary.txt")
data_summary <- read_tsv("./data/030.2_invivo_data_summary.txt")

######################################
x_labels <- c(1,2,5,10)
names(x_labels) <- c("4.3hpf","9hpf","30hpf","72hpf")
data_summary <- group_by(data_summary,dev_stage,cas9_RNP) %>%
  summarise(shannon_entopy_mean=mean(shannon_entropy),
            shannon_entropy_sd=sd(shannon_entropy),
            gini_index_mean=mean(gini_index),
            gini_index_sd=sd(gini_index),
            diversity_mean=mean(diversity),
            diversity_sd=sd(diversity),
            unmutated_percentage_mean=mean(unmutated_percentage),
            unmutated_percentage_sd=sd(unmutated_percentage)) %>%
  mutate(stage=x_labels[dev_stage]-0.2*(cas9_RNP=="0.3x"))

ggplot(data_summary,aes(x=stage,y=shannon_entopy_mean,
              ymin=shannon_entopy_mean-shannon_entropy_sd,
              ymax=shannon_entopy_mean+shannon_entropy_sd,
              color=cas9_RNP)) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  geom_point(size=5) +
  geom_errorbar(width=.02) +
  theme_classic()
ggsave("./figures/invivo_shannon_entropy.png",width=5,height=5)

ggplot(data_summary,aes(x=stage,y=gini_index_mean,
              ymin=gini_index_mean-gini_index_sd,
              ymax=gini_index_mean+gini_index_sd,
              color=cas9_RNP)) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  geom_point(size=5) +
  geom_errorbar(width=.02) +
  theme_classic()
ggsave("./figures/invivo_gini_index.png",width=5,height=5)

ggplot(data_summary,aes(x=stage,y=diversity_mean,
              ymin=diversity_mean-diversity_sd,
              ymax=diversity_mean+diversity_sd,
              color=cas9_RNP)) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  geom_point(size=5) +
  geom_errorbar(width=.02) +
  theme_classic()
ggsave("./figures/invivo_diversity.png",width=5,height=5)

ggplot(data_summary,aes(x=stage,y=unmutated_percentage_mean,
              ymin=unmutated_percentage_mean-unmutated_percentage_sd,
              ymax=unmutated_percentage_mean+unmutated_percentage_sd,
              color=cas9_RNP)) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  geom_point(size=5) +
  geom_errorbar(width=.02) +
  theme_classic()
ggsave("./figures/invivo_unmutated_percentage.png",width=5,height=5)

