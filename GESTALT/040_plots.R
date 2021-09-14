library(readr)
library(dplyr)
library(tidyr)
library(TraceQC)
library(ggplot2)
library(stringr)
library(DescTools)

ref <- parse_ref_file("./data/000_ref_v6.txt")
ref$regions <- mutate(ref$regions,region=case_when(region=="sgRNA" ~ "spacer",
                                                   TRUE ~ region))
plot_construct(ref,chr_per_row=100,chr_size=1.5)
ggsave("./figures/figs1_construct.png",width=6,height=1.5)


alignment_permutation <- read_tsv("./data/020.2_alignment_threshold_r1.txt")
plot_alignment_permutation(alignment_permutation)
ggsave("./figures/figs1_alignment_permutation.png",width=4,height=4)


model <- loess(score~permutate_percent,data=alignment_permutation)
aligned_reads <- read_tsv("./data/020.1_alignment/SRR3561150_1.txt")
ggplot(aligned_reads) +
  geom_histogram(aes(x=score,y=stat(count)/sum(count)),bins=50) +
  geom_vline(xintercept=predict(model,0.4),linetype="dashed",color="red",size=2) +
  xlab("alignment score") +
  ylab("frequency") +
  theme_classic()
ggsave("./figures/figs1_alignment_r1_score_hist.png",width=6,height=4)

pos_example_seq <- filter(aligned_reads,score>120,score<121)
pos <- data.frame(seq=strsplit(pos_example_seq$target_seq[1],split="")[[1]],
                  ref=strsplit(pos_example_seq$target_ref[1],split="")[[1]]) %>%
  mutate(x=1:n()) %>%
  mutate(is_match=(seq==ref))

neg_example_seq <- filter(aligned_reads,score>15,score<16)
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
ggsave("./figures/figs1_r1_example_sequences.png",width=12,height=2)


plot_lorenz_curve(aligned_reads)
ggsave("./figures/figs1_r1_lorenz_curve.png",width=4,height=4)


alignment_permutation <- read_tsv("./data/020.2_alignment_threshold_r2.txt")
model <- loess(score~permutate_percent,data=alignment_permutation)
aligned_reads <- read_tsv("./data/020.1_alignment/SRR3561150_2.txt")
ggplot(aligned_reads) +
  geom_histogram(aes(x=score,y=stat(count)/sum(count)),bins=50) +
  geom_vline(xintercept=predict(model,0.4),linetype="dashed",color="red",size=2) +
  xlab("alignment score") +
  ylab("frequency") +
  theme_classic()
ggsave("./figures/figs1_alignment_r2_score_hist.png",width=6,height=4)

pos_example_seq <- filter(aligned_reads,score>140,score<141)
pos <- data.frame(seq=strsplit(pos_example_seq$target_seq[1],split="")[[1]],
                  ref=strsplit(pos_example_seq$target_ref[1],split="")[[1]]) %>%
  mutate(x=1:n()) %>%
  mutate(is_match=(seq==ref))

neg_example_seq <- filter(aligned_reads,score>15,score<16)
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
ggsave("./figures/figs1_r2_example_sequences.png",width=12,height=2)


plot_lorenz_curve(aligned_reads)
ggsave("./figures/figs1_r2_lorenz_curve.png",width=4,height=4)


all_data <- readRDS("./data/030.1_all_data_whole_organism.rds")
mutations <- all_data$SRR3561150 %>%
  group_by(type,start,length,mutate_to,read) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=mean(count)) %>%
  ungroup

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


data_summary <- read_tsv("./data/040_data_summary.txt")
Chronic_induction <- filter(data_summary,str_detect(sample_name,"ChronicInduction"))
tmp <- strsplit(Chronic_induction$sample_name,"_")
treatment <- lapply(tmp,function(x) {paste(x[1:(length(x)-1)],collapse="_")}) %>%
  unlist()
time <- lapply(tmp,function(x) {x[length(x)]}) %>%
  unlist()
Chronic_induction$treatment <- treatment
Chronic_induction$time <- time
ChronicInduction_PosDox_0h <- filter(Chronic_induction,sample_name=="ChronicInduction_PosDox_0h")
Chronic_induction <- filter(Chronic_induction,treatment!="ChronicInduction_PosDox") %>%
  rbind(mutate(ChronicInduction_PosDox_0h,treatment="ChronicInduction_PosDox_High")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,treatment="ChronicInduction_PosDox_Low")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,treatment="ChronicInduction_PosDox_Medium"))

x_labels <- c(0,1,2,4,6,8)
names(x_labels)=c("0h","12h","24h","48h","72h","96h")
Chronic_induction <- mutate(Chronic_induction,group=treatment,time_int=x_labels[time])
ggplot(data=Chronic_induction,
       aes(x=time_int,y=unmutated_percentage,color=treatment)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=3,byrow=TRUE))
ggsave("./figures/figs2_ChronicInduction_unmutated.png",height=6,width=6)

ggplot(data=Chronic_induction,
       aes(x=time_int,y=shannon_entropy,color=treatment)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=3,byrow=TRUE))
ggsave("./figures/figs2_ChronicInduction_shannon_entropy.png",height=6,width=6)

