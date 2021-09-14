library(ggtern)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggalt)
library(RColorBrewer)

df <- data.frame()
################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
data_points <- data.frame()
for (sample_name in names(Carlin)[2:length(Carlin)]) {
  tmp <- filter(Carlin[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution=mutation)
df <- rbind(df,mutate(data_points,platform="Carlin"))

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(scCarlin)[2:length(scCarlin)]) {
  tmp <- filter(scCarlin[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution)
df <- rbind(df,mutate(data_points,platform="scCarlin"))

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
data_points <- data.frame()
for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  tmp <- filter(GESTALT[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution=mutation)
df <- rbind(df,mutate(data_points,platform="GESTALT"))

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(scGESTALT)[2:length(scGESTALT)]) {
  tmp <- filter(scGESTALT[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution)
df <- rbind(df,mutate(data_points,platform="scGESTALT"))

################ hgRNA-invitro ################
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- filter(hgRNA_invitro[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution=mutation)
df <- rbind(df,mutate(data_points,platform="hgRNA-invitro"))

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(hgRNA_invivo)[2:length(hgRNA_invivo)]) {
  if (nrow(hgRNA_invivo[[sample_name]])>100) {
  tmp <- filter(hgRNA_invivo[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}}
data_points <- select(data_points,deletion,insertion,substitution=mutation)
df <- rbind(df,mutate(data_points,platform="hgRNA-invivo"))

################ LINNAEUS ################
LINNAEUS <- readRDS("../LINNAEUS/data/060.1_all_data.rds")
data_points <- data.frame()
for (sample_name in names(LINNAEUS)[2:length(LINNAEUS)]) {
  tmp <- filter(LINNAEUS[[sample_name]],type!="unmutated") %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise() %>%
    ungroup %>%
    group_by(type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type) %>%
    mutate(count=count/sum(count))
  data_points <- rbind(data_points,spread(tmp,type,count))}
data_points <- select(data_points,deletion,insertion,substitution=mutation)
df <- rbind(df,mutate(data_points,platform="LINNAEUS"))


cols <- brewer.pal(7,"Set1")
names(cols) <- c("Carlin","scCarlin","GESTALT","scGESTALT",
                 "hgRNA-invitro","hgRNA-invivo","LINNAEUS")
ggtern(df,aes(deletion,insertion,substitution,fill=platform,color=platform)) +
  geom_encircle(alpha=0.2,size=1,expand=0.02) +
  geom_point(size=1,alpha=0.5) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  theme(legend.position="none",
        legend.text = element_text(size=12),
        panel.grid = element_line(color="grey"),
        panel.background = element_blank(),
        plot.background=element_blank()) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))
ggsave("./figures/fig2_ternary.png",height=4,width=4)
