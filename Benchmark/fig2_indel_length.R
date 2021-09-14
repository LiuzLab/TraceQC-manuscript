library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)

################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
data_points <- data.frame()
for (sample_name in names(Carlin)[2:length(Carlin)]) {
  tmp <- filter(Carlin[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)
labels <- c("1","2","4","8","16","32","64","128","256","256+")
# data_points$length <- factor(data_points$length,levels=factor_lvs)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_carlin.png",width=3,height=3)

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(scCarlin)[2:length(scCarlin)]) {
  tmp <- filter(scCarlin[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=n()) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_sccarlin.png",width=3,height=3)

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
data_points <- data.frame()
for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  tmp <- filter(GESTALT[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_GESTALT.png",width=3,height=3)

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(scGESTALT)[2:length(scGESTALT)]) {
  tmp <- filter(scGESTALT[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=n()) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_scGESTALT.png",width=3,height=3)

################ hgRNA-invitro ################
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- filter(hgRNA_invitro[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_hgRNA_invitro.png",width=3,height=3)

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(hgRNA_invivo)[2:length(hgRNA_invivo)]) {
  tmp <- filter(hgRNA_invivo[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("./figures/fig2_indel_length_hgRNA_invivo.png",width=3,height=3)

################ LINNAEUS ################
LINNAEUS <- readRDS("../LINNAEUS/data/060.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(LINNAEUS)[2:length(LINNAEUS)]) {
  tmp <- filter(LINNAEUS[[sample_name]],type %in% c("insertion","deletion")) %>%
    group_by(type,length) %>%
    mutate(length=ceiling(log2(length))) %>%
    mutate(case_when(length<=8 ~ length,
                     TRUE ~ 9)) %>%
    summarise(count=n()) %>%
    ungroup %>%
    group_by(type) %>%
    mutate(count=count/sum(count)) %>%
    ungroup
  data_points <- rbind(data_points,tmp)}
data_points <- group_by(data_points,type,length) %>%
  summarise(frequency=mean(count),sd=sd(count)) %>%
  ungroup %>%
  mutate(sd=replace_na(sd,0)) %>%
  mutate(is_deletion=2*(type=="deletion")-1) %>%
  mutate(frequency=frequency*is_deletion,sd=sd*is_deletion)

ggplot(data_points,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=frequency,ymax=frequency+sd),width=.2) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_x_continuous(breaks=seq(0.5,
                                max(data_points$length)+0.5),
                     labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_blank())
ggsave("./figures/fig2_indel_length_LINNAEUS.png",width=5,height=3)
