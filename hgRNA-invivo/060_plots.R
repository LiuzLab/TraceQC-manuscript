library(readr)
library(dplyr)
library(tidyr)
library(TraceQC)
library(ggplot2)
library(stringr)
library(tools)
library(DescTools)


ref <- parse_ref_file("./data/000_ref/L21.txt")
plot_construct(ref,chr_per_row=100,chr_size=1.5)
ggsave("./figures/figs1_construct.png",width=6,height=1.5)


alignment_permutation <- read_tsv("./data/030.3_alignment_threshold/alignment_threshold_L21_AAACCCCGGG.txt")
plot_alignment_permutation(alignment_permutation)
ggsave("./figures/figs1_alignment_permutation.png",width=4,height=4)


identifiers <- read_csv("./data/000_ref/hgRNA_identifiers.csv")
L21_identifiers <- filter(identifiers,Length==21) %>%
  pull(`Identifier (ID)`)

model <- loess(score~permutate_percent,data=alignment_permutation)
aligned_reads <- data.frame()
for (idf in L21_identifiers) {
  aligned_reads <- rbind(aligned_reads,
                         read_tsv(sprintf("./data/030.1_alignment/SRR7633623_%s.txt",idf)))
}

ggplot(aligned_reads) +
  geom_histogram(aes(x=score,y=stat(count)/sum(count)),bins=50) +
  geom_vline(xintercept=predict(model,0.4),linetype="dashed",color="red",size=2) +
  xlab("alignment score") +
  ylab("frequency") +
  theme_classic()
ggsave("./figures/figs1_alignment_score_hist.png",width=6,height=4)


pos_example_seq <- filter(aligned_reads,score>100,score<200)
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


all_data <- readRDS("./data/040.2_all_data_normed.rds")
mutations <- all_data$SRR7633623  %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
  
mutation_type_donut(mutations)
ggsave("./figures/figs2_donut.png",height=4,width=4)


ref <- parse_ref_file("./data/000_ref/L21_AAGCCGCGCG.txt")
mutations <- filter(mutations,identifier=="AAGCCGCGCG")
png("./figures/figs2_deletion_circular.png")
plot_deletion_hotspot(mutations,ref)
dev.off()

png("./figures/figs2_insertion_circular.png")
plot_insertion_hotspot(mutations,ref)
dev.off()

png("./figures/figs2_substitution_circular.png")
plot_point_substitution_hotspot(mutations,ref)
dev.off()


identifiers <- select(identifiers,identifier=`Identifier (ID)`,Length)
data_summary <- read_tsv("./data/050_data_summary.txt") %>%
  mutate(unmutated_count=replace_na(unmutated_count,0)) %>%
  filter((mutated_count+unmutated_count) > 800) %>%
  mutate(mutated_percentage=mutated_count/(mutated_count+unmutated_count))
  
x_labels <- c(0,3.5,6.5,7.5,8.5,10.5,12.5,16.5,22)
names(x_labels)=c("E0","E03p5","E06p5","E07p5","E08p5",
                  "E10p5","E12p5","E16p5","Adult")
df <- filter(data_summary,dev_stage!="P21") %>%
  group_by(dev_stage,identifier) %>%
  summarise(shannon_entopy_mean=mean(SE),
            shannon_entropy_sd=sd(SE),
            mutated_percentage_mean=mean(mutated_percentage),
            mutated_percentage_sd=sd(mutated_percentage),
  ) %>%
  mutate(stage=x_labels[dev_stage])

df <- rbind(df,data.frame(dev_stage="E0",identifier=unique(df$identifier),
           shannon_entopy_mean=0,
           shannon_entropy_sd=0,
           mutated_percentage_mean=0,
           mutated_percentage_sd=0,
           stage=0)) %>%
  left_join(identifiers) %>%
  mutate(Length=as.factor(Length))

names(x_labels)=c("E0","E3.5","E6.5","E7.5","E8.5",
                  "E10.5","E12.5","E16.5","Adult")
ggplot(df,aes(x=stage,y=shannon_entopy_mean,group=identifier,
              ymin=shannon_entopy_mean-shannon_entropy_sd,
              ymax=shannon_entopy_mean+shannon_entropy_sd)) +
  geom_line(color='#E69F00') +
  scale_x_continuous(breaks=x_labels) +
  geom_point(size=5,color='#E69F00') +
  geom_errorbar(width=.02,color='#E69F00') +
  theme_classic()
ggsave("./figures/figs3_shannon_entropy.png",width=5,height=5)

ggplot(filter(df,Length==21),aes(x=stage,y=mutated_percentage_mean,group=identifier)) +
  geom_line(alpha=0.5,color='#E69F00') +
  scale_x_continuous(breaks=x_labels) +
  geom_point(alpha=0.5,size=3,color='#E69F00') +
  ggtitle("L21") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=30,vjust=0.5,size=8))
ggsave("./figures/figs3_L21_mutated_percentage.png",width=4,height=4)

ggplot(filter(df,Length==25),aes(x=stage,y=mutated_percentage_mean,group=identifier)) +
  geom_line(alpha=0.5,color='#E69F00') +
  scale_x_continuous(breaks=x_labels) +
  geom_point(alpha=0.5,size=3,color='#E69F00') +
  ggtitle("L25") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=30,vjust=0.5,size=8))
ggsave("./figures/figs3_L25_mutated_percentage.png",width=4,height=4)

ggplot(filter(df,Length==30),aes(x=stage,y=mutated_percentage_mean,group=identifier)) +
  geom_line(alpha=0.5,color='#E69F00') +
  scale_x_continuous(breaks=x_labels) +
  geom_point(alpha=0.5,size=3,color='#E69F00') +
  ggtitle("L30") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=30,vjust=0.5,size=8))
ggsave("./figures/figs3_L30_mutated_percentage.png",width=4,height=4)

ggplot(filter(df,Length==35),aes(x=stage,y=mutated_percentage_mean,group=identifier)) +
  geom_line(alpha=0.5,color='#E69F00') +
  scale_x_continuous(breaks=x_labels) +
  geom_point(alpha=0.5,size=3,color='#E69F00') +
  ggtitle("L35") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=30,vjust=0.5,size=8))
ggsave("./figures/figs3_L35_mutated_percentage.png",width=4,height=4)

