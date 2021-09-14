library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)

hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")

df1 <- hgRNA_invitro$`A21-14d(A')` %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep1_count=sum(count)) %>%
  ungroup

df2 <- hgRNA_invitro$`A21-5d` %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep2_count=sum(count)) %>%
  ungroup

df3 <- hgRNA_invitro$`A21-pop6` %>%
  filter(type!="unmutated") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep3_count=sum(count)) %>%
  ungroup

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep1_count,tmp$rep2_count,method="pearson"),2)
  return(c(per,pearson))}

df <- full_join(df1,df2,by=c("type","start","length","mutate_to")) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep2_count=log10(rep2_count+1)) %>%
  mutate(missing=(rep1_count==0|rep2_count==0)) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

ggplot(df,aes(x=rep1_count,y=rep2_count,color=type,shape=type)) +
  geom_point(size=0.8,alpha=0.4) +
  coord_cartesian() +
  scale_color_manual(values=c("deletion"="#F8766D",
                             "insertion"="#00BFC4",
                             "substitution"="#7CAE00")) +
  annotate("text",x=4,y=1,size=2,
           label=sprintf("deletion %s %s\n insertion %s %s\n substitution %s %s\n",
                    scores(filter(df,type=="deletion"))[1],
                    scores(filter(df,type=="deletion"))[2],
                    scores(filter(df,type=="insertion"))[1],
                    scores(filter(df,type=="insertion"))[2],
                    scores(filter(df,type=="substitution"))[1],
                    scores(filter(df,type=="substitution"))[2])) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    legend.position = "none")
ggsave("./figures/figs12_cor_rep12_hgRNA_invitro.png",width=3,height=3)

df <- full_join(df1,df3,by=c("type","start","length","mutate_to")) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep1_count==0|rep3_count==0)) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep1_count,tmp$rep3_count,method="pearson"),2)
  return(c(per,pearson))}

ggplot(df,aes(x=rep1_count,y=rep3_count,color=type,shape=type)) +
  geom_point(size=0.8,alpha=0.4) +
  coord_cartesian() +
  scale_color_manual(values=c("deletion"="#F8766D",
                              "insertion"="#00BFC4",
                              "substitution"="#7CAE00")) +
  annotate("text",x=4,y=1,size=2,
           label=sprintf("deletion %s %s\n insertion %s %s\n substitution %s %s\n",
                         scores(filter(df,type=="deletion"))[1],
                         scores(filter(df,type=="deletion"))[2],
                         scores(filter(df,type=="insertion"))[1],
                         scores(filter(df,type=="insertion"))[2],
                         scores(filter(df,type=="substitution"))[1],
                         scores(filter(df,type=="substitution"))[2])) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    legend.position = "none")
ggsave("./figures/figs12_cor_rep13_hgRNA_invitro.png",width=3,height=3)


df <- full_join(df2,df3,by=c("type","start","length","mutate_to")) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep2_count=log10(rep2_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep2_count==0|rep3_count==0)) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

scores <- function (df) {
  per <- round(1-sum(df$missing)/nrow(df),2)
  tmp <- filter(df,!missing)
  pearson <- round(cor(tmp$rep2_count,tmp$rep3_count,method="pearson"),2)
  return(c(per,pearson))}

ggplot(df,aes(x=rep2_count,y=rep3_count,color=type,shape=type)) +
  geom_point(size=0.8,alpha=0.4) +
  coord_cartesian() +
  scale_color_manual(values=c("deletion"="#F8766D",
                              "insertion"="#00BFC4",
                              "substitution"="#7CAE00")) +
  annotate("text",x=4,y=1,size=2,
           label=sprintf("deletion %s %s\n insertion %s %s\n substitution %s %s\n",
                         scores(filter(df,type=="deletion"))[1],
                         scores(filter(df,type=="deletion"))[2],
                         scores(filter(df,type=="insertion"))[1],
                         scores(filter(df,type=="insertion"))[2],
                         scores(filter(df,type=="substitution"))[1],
                         scores(filter(df,type=="substitution"))[2])) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    legend.position = "none")
ggsave("./figures/figs12_cor_rep23_hgRNA_invitro.png",width=3,height=3)


data_summary <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- strsplit(sample_name,"-") %>% unlist()
  barcode <- tmp[1]
  day <- tmp[2]
  if (str_detect(sample_name,"\\(A'\\)")) {
    barcode <- "A'21"
    day <- substr(day,start=1,stop=nchar(day)-4)}
  df <- hgRNA_invitro[[sample_name]] %>%
    filter(type!="unmutated") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    group_by(type) %>%
    summarise(amount=n(),shannon_entropy=Entropy(count)) %>%
    ungroup
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name,
                               barcode=barcode,time=day,))}

x_labels <- c(0,2,14)
names(x_labels)=c("0d","2d","14d")
data_summary <- filter(data_summary,barcode=="A\'21") %>%
  mutate(time=x_labels[time]) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

ggplot(data=data_summary,aes(x=time,y=amount,shape=type,
                   group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("deletion"="#F8766D",
                              "insertion"="#00BFC4",
                              "substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"2d"=2,"14d"=14)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    legend.position = "none")
ggsave("./figures/figs12_lineplot_amount_hgRNA_invitro.png",height=3,width=3)

ggplot(data=data_summary,aes(x=time,y=shannon_entropy,shape=type,
                             group=type,color=type)) +
  geom_point(size=4,alpha=0.8) +
  geom_line() +
  scale_color_manual(values=c("deletion"="#F8766D",
                              "insertion"="#00BFC4",
                              "substitution"="#7CAE00")) +
  scale_x_continuous(breaks=c("0"=0,"2d"=2,"14d"=14)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.bottom = element_text(size=15),
    legend.position = "none")
ggsave("./figures/figs12_lineplot_shannon_entropy_hgRNA_invitro.png",height=3,width=3)
