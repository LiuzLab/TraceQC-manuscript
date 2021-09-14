library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)

############## Carlin ##############
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
data_summary <- data.frame()
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2

for (sample_name in names(Carlin)[2:length(Carlin)]) {
  df <- Carlin[[sample_name]] %>%
    filter(type!="unmutated") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
           ending_sgRNA=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                          (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                          TRUE ~ type)) %>%
    mutate(starting_cutsite=start+1-target_end[starting_sgRNA],
           ending_cutsite=start+length+1-target_end[ending_sgRNA]) %>%
    group_by(type) %>%
    summarise(amount=n(),shannon_entropy=Entropy(count)) %>%
    ungroup
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}

Chronic_induction <- filter(data_summary,str_detect(sample_name,"ChronicInduction_PosDox"))
tmp <- strsplit(Chronic_induction$sample_name,"_")
x_labels <- c(0,3,6,12,18,24)
names(x_labels)=c("0h","12h","24h","48h","72h","96h")
group <- lapply(tmp,function(x) {paste(x[1:(length(x)-1)],collapse="_")}) %>%
  unlist()
time <- lapply(tmp,function(x) {x[length(x)]}) %>%
  unlist()
Chronic_induction$group <- group
Chronic_induction$time <- x_labels[time]

ChronicInduction_PosDox_0h <- filter(Chronic_induction,sample_name=="ChronicInduction_PosDox_0h")
Chronic_induction <- filter(Chronic_induction,group!="ChronicInduction_PosDox") %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_High")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_Low")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_Medium")) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

df <- filter(Chronic_induction,group=="ChronicInduction_PosDox_Low")
ggplot(data=df,aes(x=time,y=amount,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                       "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig5_lineplot_amount_carlin_low_DOX.png",height=3,width=3)

ggplot(data=df,aes(x=time,y=shannon_entropy,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig5_lineplot_SE_carlin_low_DOX.png",height=3,width=3)

df <- filter(Chronic_induction,group=="ChronicInduction_PosDox_Medium")
ggplot(data=df,aes(x=time,y=amount,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/figs11_lineplot_amount_carlin_medium_DOX.png",height=3,width=3)

ggplot(data=df,aes(x=time,y=shannon_entropy,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/figs11_lineplot_SE_carlin_medium_DOX.png",height=3,width=3)

df <- filter(Chronic_induction,group=="ChronicInduction_PosDox_High")
ggplot(data=df,aes(x=time,y=amount,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/figs11_lineplot_amount_carlin_high_DOX.png",height=3,width=3)

ggplot(data=df,aes(x=time,y=shannon_entropy,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/figs11_lineplot_SE_carlin_high_DOX.png",height=3,width=3)
