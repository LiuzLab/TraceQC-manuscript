library(readr)
library(dplyr)
library(tidyr)
library(TraceQC)
library(ggplot2)
library(stringr)
library(DescTools)

ref <- parse_ref_file("./data/000_ref.txt")
target_start <- filter(ref$regions,region=="target") %>% pull(start)
target_end <- filter(ref$regions,region=="target") %>% pull(end)
sgRNA_start <- filter(ref$regions,region=="sgRNA") %>% pull(start)
sgRNA_end <- filter(ref$regions,region=="sgRNA") %>% pull(end)
ref$regions <- mutate(ref$regions,region=case_when(region=="sgRNA" ~ "spacer",
                                                   TRUE ~ region))
ref$regions$start <- c(sgRNA_start-20,sgRNA_start)

plot_construct(ref,chr_per_row=100,chr_size=1.5)
ggsave("./figures/figs1_construct.png",width=6,height=1.5)

alignment_permutation <- read_tsv("./data/031_alignment_threshold.txt")
plot_alignment_permutation(alignment_permutation)
ggsave("./figures/figs1_alignment_permutation.png",width=4,height=4)

model <- loess(score~permutate_percent,data=alignment_permutation)
aligned_reads <- read_tsv("./data/030_traceQC_alignment/SRR6211489.txt")
ggplot(aligned_reads) +
  geom_histogram(aes(x=score,y=stat(count)),bins=50) +
  geom_vline(xintercept=predict(model,0.4),linetype="dashed",color="red",size=2) +
  scale_y_continuous(trans='log10') +
  xlab("alignment score") +
  ylab("count") +
  theme_classic()
ggsave("./figures/figs1_alignment_score_hist.png",width=6,height=4)


pos_example_seq <- filter(aligned_reads,score>59,score<60)
pos <- data.frame(seq=strsplit(pos_example_seq$target_seq[1],split="")[[1]],
                  ref=strsplit(pos_example_seq$target_ref[1],split="")[[1]]) %>%
  mutate(x=1:n()) %>%
  mutate(is_match=(seq==ref))

neg_example_seq <- filter(aligned_reads,score>20,score<21)
neg <- data.frame(seq=strsplit(neg_example_seq$target_seq[1],split="")[[1]],
                  ref=strsplit(neg_example_seq$target_ref[1],split="")[[1]]) %>%
  mutate(x=1:n()) %>%
  mutate(is_match=(seq==ref))
ggplot(data.frame()) +
  geom_text(aes(x=x,y=6,label=seq,color=is_match),data=pos,size=1.5) +
  geom_text(aes(x=x,y=5,label=ref,color=is_match),data=pos,size=1.5) +
  geom_text(aes(x=x,y=1,label=seq,color=is_match),data=neg,size=1.5) +
  geom_text(aes(x=x,y=0,label=ref,color=is_match),data=neg,size=1.5) +
  scale_color_manual(values=c("TRUE"="black","FALSE"="red")) +
  theme_void()
ggsave("./figures/figs1_example sequences.png",width=12,height=2)


aligned_reads <- filter(aligned_reads,score>predict(model,0.4))
plot_lorenz_curve(aligned_reads)
ggsave("./figures/figs1_lorenz_curve.png",width=4,height=4)


UMI_count <- get_UMI_count_per_CB(aligned_reads) %>%
  filter(CB != "/") %>%
  arrange(desc(UMI_per_CB)) %>%
  mutate(x=1:n())
ggplot(UMI_count,aes(x=x,y=UMI_per_CB)) + 
  geom_line(color="red",size=2) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab("Cells") +
  ylab("UMI count") +
  # scale_y_continuous(breaks=seq(0,1,0.25)) +
  theme(axis.text.x = element_text(size=15),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(size=15))
ggsave("./figures/fig1_umi_count.png",width=4,height=4)

read_count_per_UMI <- get_read_count_per_UB(aligned_reads) %>%
  filter(UB != "/") %>%
  filter(CB != "/") %>%
  arrange(desc(read_count_per_UMI)) %>%
  mutate(x=1:n())
ggplot(read_count_per_UMI,aes(x=x,y=read_count_per_UMI)) + 
  geom_line(color="red",size=2) +
  # coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab("UMI") +
  ylab("Read count") +
  # scale_y_continuous(breaks=seq(0,1,0.25)) +
  theme(axis.text.x = element_text(size=15),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(size=15))
ggsave("./figures/fig1_umi_read_count.png",width=4,height=4)


all_data <- readRDS("./data/060.1_all_data.rds")
mutations <- all_data$SRR6211489
mutation_type_donut(mutations)
ggsave("./figures/figs2_donut.png",height=4,width=4)

mutations <- group_by(mutations,type,start,length,mutate_to) %>%
  summarise(count=n()) %>%
  ungroup


png("./figures/figs2_deletion_circular.png")
plot_deletion_hotspot(mutations,ref,use_log_count=FALSE)
dev.off()

png("./figures/figs2_insertion_circular.png")
plot_insertion_hotspot(mutations,ref,use_log_count=FALSE)
dev.off()

png("./figures/figs2_substitution_circular.png")
plot_point_substitution_hotspot(mutations,ref)
dev.off()


df <- data.frame()
for (sample_name in names(all_data)[2:length(names(all_data))]) {
  total_numer_of_cells <- length(unique(all_data[[sample_name]]$CB))
  tmp <- all_data[[sample_name]] %>%
    filter(type!="unmutated") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(number_of_cells = n()) %>%
    ungroup %>%
    mutate(percentage=number_of_cells/total_numer_of_cells) %>%
    arrange(desc(percentage)) %>%
    mutate(barcode_rank=1:n()) %>%
    mutate(sample_name=sample_name)
  df <- rbind(df,tmp)}

ggplot(df,aes(x=barcode_rank,y=number_of_cells,group=sample_name)) + 
  geom_line(color="red",alpha=0.5) +
  # coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab("Mutations") +
  ylab("Number of cells") +
  # scale_y_continuous(breaks=seq(0,1,0.25)) +
  theme(axis.text.x = element_text(size=15),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(size=15))
ggsave("./figures/figs2_mutation_distribution.png",height=4,width=4)
