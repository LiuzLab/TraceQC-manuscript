library(readr)
library(dplyr)
library(ggplot2)

Carlin <- readRDS("../Carlin/data/030_all_data.rds")
df1 <- Carlin$ChronicInduction_PosDox_Low_96h %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep1_count=sum(count)) %>%
  ungroup

df2 <- Carlin$ChronicInduction_PosDox_Medium_96h %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep2_count=sum(count)) %>%
  ungroup

df3 <- Carlin$ChronicInduction_PosDox_High_96h %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep3_count=sum(count)) %>%
  ungroup

plotting_df <- data.frame()
df <- full_join(df1,df2) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep2_count=log10(rep2_count+1)) %>%
  mutate(missing=(rep1_count==0|rep2_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                  TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep12")
plotting_df <- rbind(plotting_df,df)

df <- full_join(df1,df3) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep1_count==0|rep3_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                             TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep13")
plotting_df <- rbind(plotting_df,df)

df <- full_join(df2,df3) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep2_count=log10(rep2_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep2_count==0|rep3_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                             TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep23")
plotting_df <- rbind(plotting_df,df)

ggplot(plotting_df,aes(x=log_length,y=overlap)) +
  geom_smooth(method="loess",se = FALSE, color = "grey", size=0.5) +
  geom_point(aes(color=rep),alpha=0.5) +
  xlab("log(insertion length)") +
  ylab("percentage of identical mutations") +
  # geom_line(alpha=0.5) +
  coord_cartesian() +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    # axis.title = element_blank(),
    legend.position = "bottom")
ggsave("./figures/figs10_insertion_length_carlin.png",width=4,height=4)


hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
df1 <- hgRNA_invitro$`A21-14d(A')` %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep1_count=sum(count)) %>%
  ungroup

df2 <- hgRNA_invitro$`A21-5d` %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep2_count=sum(count)) %>%
  ungroup

df3 <- hgRNA_invitro$`A21-pop6` %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(rep3_count=sum(count)) %>%
  ungroup

plotting_df <- data.frame()
df <- full_join(df1,df2) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep2_count=log10(rep2_count+1)) %>%
  mutate(missing=(rep1_count==0|rep2_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                  TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep12")
plotting_df <- rbind(plotting_df,df)

df <- full_join(df1,df3) %>%
  mutate(rep1_count=replace_na(rep1_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep1_count=log10(rep1_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep1_count==0|rep3_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                             TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep13")
plotting_df <- rbind(plotting_df,df)

df <- full_join(df2,df3) %>%
  mutate(rep2_count=replace_na(rep2_count,0)) %>%
  mutate(rep3_count=replace_na(rep3_count,0)) %>%
  mutate(rep2_count=log10(rep2_count+1),
         rep3_count=log10(rep3_count+1)) %>%
  mutate(missing=(rep2_count==0|rep3_count==0)) %>%
  # mutate(log_length=ceiling(log2(length))) %>%
  mutate(log_length=log2(length)) %>%
  # mutate(log_length=case_when(log_length<=8 ~ log_length,
  #                             TRUE ~ 9)) %>%
  group_by(log_length) %>%
  summarise(overlap=sum(!missing)/n()) %>%
  ungroup %>%
  mutate(rep="rep23")
plotting_df <- rbind(plotting_df,df)

ggplot(plotting_df,aes(x=log_length,y=overlap)) +
  geom_smooth(method="loess",se = FALSE, color = "grey", size=0.5) +
  geom_point(aes(color=rep),alpha=0.5) +
  xlab("log(insertion length)") +
  ylab("percentage of identical mutations") +
  # geom_line(alpha=0.5) +
  coord_cartesian() +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    # axis.title = element_blank(),
    legend.position = "bottom")
ggsave("./figures/figs10_insertion_length_hgRNA_invitro.png",width=4,height=4)
