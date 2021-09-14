library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)

############## Carlin ##############
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- ref_link_end-ref_target_start+2

df1 <- Carlin$ChronicInduction_NegDox_Trial1_0h %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length) %>%
  summarise(rep1_count=sum(count)) %>%
  ungroup

df2 <- Carlin$ChronicInduction_NegDox_Trial2_0h %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length) %>%
  summarise(rep2_count=sum(count)) %>%
  ungroup

df3 <- Carlin$ChronicInduction_PosDox_0h %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length) %>%
  summarise(rep3_count=sum(count)) %>%
  ungroup

df <- full_join(df1,df2,by=c("type","start","length")) %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end),
         ending_sgRNA=findInterval(start+length,ref_link_end)) %>%
  mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                        (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                        .data$type=="mutation" ~ "substitution",
                        TRUE ~ type)) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep2_count=log10(rep2_count+1)) %>%
  mutate(missing=(rep1_count==0|rep2_count==0))

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep1_count,tmp$rep2_count,method="pearson"),2)
  return(c(per,pearson))}
scores(filter(df,type=="inter-target deletion"))
scores(filter(df,type=="intra-target deletion"))
scores(filter(df,type=="insertion"))
scores(filter(df,type=="substitution"))

ggplot(df,aes(x=rep1_count,y=rep2_count,color=type,shape=type)) +
  geom_point(size=3,alpha=0.6) +
  coord_cartesian() +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig5_cor_rep12_carlin.png",width=3,height=3)


df <- full_join(df1,df3,by=c("type","start","length")) %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end),
         ending_sgRNA=findInterval(start+length,ref_link_end)) %>%
  mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                        (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                        .data$type=="mutation" ~ "substitution",
                        TRUE ~ type)) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep1_count==0|rep3_count==0))

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep1_count,tmp$rep3_count,method="pearson"),2)
  return(c(per,pearson))}

ggplot(df,aes(x=rep1_count,y=rep3_count,color=type,shape=type)) +
  geom_point(size=3,alpha=0.6) +
  coord_cartesian() +
  annotate("text",x=4,y=1,size=2,
           label=sprintf("inter deletion %s %s\n intra deletion %s %s\n insertion %s %s\n substitution %s %s\n",
                         scores(filter(df,type=="inter-target deletion"))[1],
                         scores(filter(df,type=="inter-target deletion"))[2],
                         scores(filter(df,type=="intra-target deletion"))[1],
                         scores(filter(df,type=="intra-target deletion"))[2],
                         scores(filter(df,type=="insertion"))[1],
                         scores(filter(df,type=="insertion"))[2],
                         scores(filter(df,type=="substitution"))[1],
                         scores(filter(df,type=="substitution"))[2])) +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig5_cor_rep13_carlin.png",width=3,height=3)


df <- full_join(df2,df3,by=c("type","start","length")) %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end),
         ending_sgRNA=findInterval(start+length,ref_link_end)) %>%
  mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                        (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                        .data$type=="mutation" ~ "substitution",
                        TRUE ~ type)) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep2_count=log10(rep2_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep2_count==0|rep3_count==0))

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep2_count,tmp$rep3_count,method="pearson"),2)
  return(c(per,pearson))}

ggplot(df,aes(x=rep2_count,y=rep3_count,color=type,shape=type)) +
  geom_point(size=3,alpha=0.6) +
  coord_cartesian() +
  annotate("text",x=4,y=1,size=2,
           label=sprintf("inter deletion %s %s\n intra deletion %s %s\n insertion %s %s\n substitution %s %s\n",
                         scores(filter(df,type=="inter-target deletion"))[1],
                         scores(filter(df,type=="inter-target deletion"))[2],
                         scores(filter(df,type=="intra-target deletion"))[1],
                         scores(filter(df,type=="intra-target deletion"))[2],
                         scores(filter(df,type=="insertion"))[1],
                         scores(filter(df,type=="insertion"))[2],
                         scores(filter(df,type=="substitution"))[1],
                         scores(filter(df,type=="substitution"))[2])) +
  scale_color_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig5_cor_rep23_carlin.png",width=3,height=3)
