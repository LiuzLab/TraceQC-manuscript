library(readr)
library(dplyr)
library(tidyr)
library(TraceQC)
library(ggplot2)
library(stringr)
library(DescTools)

ref <- parse_ref_file("./data/000_ref/A21.txt")
ref$regions <- mutate(ref$regions,region=case_when(region=="sgRNA" ~ "spacer",
                                                   TRUE ~ region))
plot_construct(ref,chr_per_row=100,chr_size=1.5)
ggsave("./figures/figs1_construct.png",width=6,height=1.5)


alignment_permutation <- read_tsv("./data/030.1_alignment_threshold/alignment_threshold_A21.txt")
plot_alignment_permutation(alignment_permutation)
ggsave("./figures/figs1_alignment_permutation.png",width=4,height=4)


model <- loess(score~permutate_percent,data=alignment_permutation)
aligned_reads <- read_tsv("./data/020_alignment/A21-14d(A').txt")
ggplot(aligned_reads) +
  geom_histogram(aes(x=score,y=stat(count)/sum(count)),bins=50) +
  geom_vline(xintercept=predict(model,0.4),linetype="dashed",color="red",size=2) +
  xlab("alignment score") +
  ylab("frequency") +
  theme_classic()
ggsave("./figures/figs1_alignment_score_hist.png",width=6,height=4)


pos_example_seq <- filter(aligned_reads,score>100,score<111)
pos <- data.frame(seq=strsplit(pos_example_seq$target_seq[1],split="")[[1]],
                  ref=strsplit(pos_example_seq$target_ref[1],split="")[[1]]) %>%
  mutate(x=1:n()) %>%
  mutate(is_match=(seq==ref))

neg_example_seq <- filter(aligned_reads,score>0,score<21)
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
ggsave("./figures/figs1_example_sequences.png",width=12,height=2)


plot_lorenz_curve(aligned_reads)
ggsave("./figures/figs1_lorenz_curve.png",width=4,height=4)


all_data <- readRDS("./data/040.1_all_data.rds")
mutations <- all_data$`A21-14d(A')` %>%
    mutate(type=case_when(type=="mutation" ~ "substitution",
                          TRUE ~ type))
mutation_type_donut(mutations)
ggsave("./figures/figs2_donut.png",height=4,width=4)


png("./figures/figs2_deletion_circular.png")
plot_deletion_hotspot(mutations,ref)
dev.off()

png("./figures/figs2_insertion_circular.png")
plot_insertion_hotspot(mutations,ref)
dev.off()

png("./figures/figs2_substitution_circular.png")
plot_point_substitution_hotspot(mutations,ref)
dev.off()


data_summary <- read_tsv("./data/040.2_data_summary.txt") %>%
  filter(!str_detect(time,"pop"))
x_labels <- c(0,1,5,1,3,5)
names(x_labels)=c("0d","2d","14d","1d","3d","5d")
data_summary <- mutate(data_summary,group=barcode,
                       plotting_group=case_when(str_detect(barcode,"A'") ~ "g1",
                                                TRUE ~ "g2")) %>%
  mutate(time_int=x_labels[time])

g1 <- filter(data_summary,plotting_group=="g1")
ggplot(data=g1,aes(x=time_int,y=unmutated_percentage,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("unmutated percentage") +
  scale_x_continuous(breaks=x_labels[1:3]) +
  theme_classic()
ggsave("./figures/figs3_g1_unmutated.png",height=4,width=4)

ggplot(data=g1,aes(x=time_int,y=shannon_entropy,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("Shannon entropy") +
  scale_x_continuous(breaks=x_labels[1:3]) +
  theme_classic()
ggsave("./figures/figs3_g1_shannon_entropy.png",height=4,width=4)

ggplot(data=g1,aes(x=time_int,y=avg_deletion_length,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("Avg Deletion Length") +
  scale_x_continuous(breaks=x_labels[1:3]) +
  theme_classic()
ggsave("./figures/figs3_g1_del_length.png",height=4,width=4)

g2 <- filter(data_summary,plotting_group=="g2")
ggplot(data=g2,aes(x=time_int,y=unmutated_percentage,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("unmutated percentage") +
  scale_x_continuous(breaks=x_labels[c(1,4:6)]) +
  theme_classic()
ggsave("./figures/figs3_g2_unmutated.png",height=4,width=4)

ggplot(data=g2,aes(x=time_int,y=shannon_entropy,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("Shannon entropy") +
  scale_x_continuous(breaks=x_labels[c(1,4:6)]) +
  theme_classic()
ggsave("./figures/figs3_g2_shannon_entropy.png",height=4,width=4)

ggplot(data=g2,aes(x=time_int,y=avg_deletion_length,
                   group=barcode,color=group)) +
  geom_point(size=5) +
  geom_line() +
  xlab("time") +
  ylab("Avg Deletion Length") +
  scale_x_continuous(breaks=x_labels[1:3]) +
  theme_classic()
ggsave("./figures/figs3_g2_del_length.png",height=4,width=4)
