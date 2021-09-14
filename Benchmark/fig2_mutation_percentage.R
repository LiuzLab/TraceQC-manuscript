library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(DescTools)
library(ggplot2)
library(TraceQC)
library(purrr)
library(tools)

df <- data.frame()
################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref <- filter(ref$regions,region=="target")
ref_length <- ref$end-ref$start+1
data_points <- data.frame()
for (sample_name in names(Carlin)[2:length(Carlin)]) {
  tmp <- filter(Carlin[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  avg_insertion_length <- sum(tmp$insertion*tmp$count)/sum(tmp$count)
  avg_substitution_length <- sum(tmp$substitution*tmp$count)/sum(tmp$count)
  data_points <- rbind(data_points,
        data.frame(deletion=avg_deletion_length/ref_length,
                   insertion=avg_insertion_length/ref_length,
                   substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="Carlin"))

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
ref <- parse_ref_file("../scCarlin/data/000_ref.txt")
ref <- filter(ref$regions,region=="target")
ref_length <- ref$end-ref$start+1
data_points <- data.frame()
max_length <- function(df) {
  res <- rep(0,ref_length)
  start <- df$start
  end <- df$start + df$length - 1
  for (i in 1:length(start)) {
    res[start[i]:end[i]] <- 1} 
  return(sum(res))}

for (sample_name in names(scCarlin)[2:length(scCarlin)]) {
  tmp <- filter(scCarlin[[sample_name]],type!="unmutated") %>%
    group_by(CB,type) %>%
    nest() %>%
    mutate(length=unlist(map(data,max_length))) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  
  avg_deletion_length <- sum(tmp$deletion)/nrow(tmp)
  avg_insertion_length <- sum(tmp$insertion)/nrow(tmp)
  avg_substitution_length <- sum(tmp$substitution)/nrow(tmp)
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length/ref_length,
                                  insertion=avg_insertion_length/ref_length,
                                  substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="scCarlin"))

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
ref <- filter(ref$regions,region=="target")
ref_length <- ref$end-ref$start+1
data_points <- data.frame()
for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  tmp <- filter(GESTALT[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  avg_insertion_length <- sum(tmp$insertion*tmp$count)/sum(tmp$count)
  avg_substitution_length <- sum(tmp$substitution*tmp$count)/sum(tmp$count)
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length/ref_length,
                                  insertion=avg_insertion_length/ref_length,
                                  substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="GESTALT"))

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
ref <- parse_ref_file("../scGESTALT/data/000_ref.txt")
ref <- filter(ref$regions,region=="target")
ref_length <- ref$end-ref$start+1
data_points <- data.frame()
max_length <- function(df) {
  res <- rep(0,ref_length)
  start <- df$start
  end <- df$start + df$length - 1
  for (i in 1:length(start)) {
    res[start[i]:end[i]] <- 1} 
  return(sum(res))}

for (sample_name in names(scGESTALT)[2:length(scGESTALT)]) {
  tmp <- filter(scGESTALT[[sample_name]],type!="unmutated") %>%
    group_by(CB,type) %>%
    nest() %>%
    mutate(length=unlist(map(data,max_length))) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  
  avg_deletion_length <- sum(tmp$deletion)/nrow(tmp)
  avg_insertion_length <- sum(tmp$insertion)/nrow(tmp)
  avg_substitution_length <- sum(tmp$substitution)/nrow(tmp)
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length/ref_length,
                                  insertion=avg_insertion_length/ref_length,
                                  substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="scGESTALT"))

################ hgRNA-invitro ################
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
all_ref_length <- list()
for (ref_file in list.files("../hgRNA-invitro/data/000_ref/")) {
  ref <- parse_ref_file(sprintf("../hgRNA-invitro/data/000_ref/%s",ref_file))
  ref <- filter(ref$regions,region=="target")
  all_ref_length[[file_path_sans_ext(ref_file)]] <- ref$end-ref$start+1}

data_points <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- filter(hgRNA_invitro[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           mutation=replace_na(mutation,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  avg_insertion_length <- sum(tmp$insertion*tmp$count)/sum(tmp$count)
  avg_substitution_length <- sum(tmp$mutation*tmp$count)/sum(tmp$count)
  ref_length <- all_ref_length[[strsplit(sample_name,split="-")[[1]][1]]]
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length/ref_length,
                                  insertion=avg_insertion_length/ref_length,
                                  substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="hgRNA_invitro"))

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
all_ref_length <- c()
identifiers <- c()
for (ref_file in list.files("../hgRNA-invivo/data/000_ref/",pattern="L[0-9][0-9]_")) {
  identifier <- file_path_sans_ext(ref_file)
  identifiers <- c(identifiers,strsplit(identifier,split="_")[[1]][2])
  ref <- parse_ref_file(sprintf("../hgRNA-invivo/data/000_ref/%s",ref_file))
  ref <- filter(ref$regions,region=="target")
  all_ref_length <- c(all_ref_length,ref$end-ref$start+1)}
names(all_ref_length) <- identifiers

data_points <- data.frame()
for (sample_name in names(hgRNA_invivo)[2:length(hgRNA_invivo)]) {
  if (nrow(hgRNA_invivo[[sample_name]])>100) {
  tmp <- filter(hgRNA_invivo[[sample_name]],type!="unmutated") %>%
    mutate(length=length/all_ref_length[identifier]) %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           mutation=replace_na(mutation,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  avg_insertion_length <- sum(tmp$insertion*tmp$count)/sum(tmp$count)
  avg_substitution_length <- sum(tmp$mutation*tmp$count)/sum(tmp$count)
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length,
                                  insertion=avg_insertion_length,
                                  substitution=avg_substitution_length))}}
df <- rbind(df,mutate(data_points,platform="hgRNA-invivo"))

################ LINNAEUS ################
LINNAEUS <- readRDS("../LINNAEUS/data/060.1_all_data.rds")
ref <- parse_ref_file("../LINNAEUS/data/000_ref.txt")
ref <- filter(ref$regions,region=="sgRNA")
ref_length <- ref$end-ref$start+41
data_points <- data.frame()
max_length <- function(df) {
  res <- rep(0,ref_length)
  start <- df$start
  end <- df$start + df$length - 1
  for (i in 1:length(start)) {
    res[start[i]:end[i]] <- 1} 
  return(sum(res))}

for (sample_name in names(LINNAEUS)[2:length(LINNAEUS)]) {
  tmp <- filter(LINNAEUS[[sample_name]],type!="unmutated") %>%
    mutate(start=start-(ref$start-20)) %>%
    group_by(CB,type) %>%
    nest() %>%
    mutate(length=unlist(map(data,max_length))) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  
  avg_deletion_length <- sum(tmp$deletion)/nrow(tmp)
  avg_insertion_length <- sum(tmp$insertion)/nrow(tmp)
  avg_substitution_length <- sum(tmp$substitution)/nrow(tmp)
  data_points <- rbind(data_points,
                       data.frame(deletion=avg_deletion_length/ref_length,
                                  insertion=avg_insertion_length/ref_length,
                                  substitution=avg_substitution_length/ref_length))}
df <- rbind(df,mutate(data_points,platform="LINNAEUS"))

################ Stacked bar plot ################
plotting_df <- gather(df,type,percentage,
                      c("deletion","insertion","substitution")) %>%
  group_by(platform,type) %>%
  summarise(percentage=mean(percentage)) %>%
  ungroup

order <- group_by(plotting_df,platform) %>%
  summarise(percentage=sum(percentage)) %>%
  ungroup %>%
  arrange(percentage) %>%
  pull(platform)
ggplot(plotting_df,aes(x=factor(platform,levels=order),y=percentage,
                       fill=factor(type,levels=c("deletion","insertion","substitution")))) + 
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1)) +
  coord_flip() +
  scale_fill_manual(values=c("deletion"="#F8766D",
                              "insertion"="#00BFC4",
                              "substitution"="#7CAE00")) +
  theme(
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color="black"),
        legend.position="bottom",
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.line.y = element_blank(),
        # axis.text.x = element_text(angle=315),
        axis.title = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_blank())
ggsave("./figures/fig2_mutation_percentage_stacked_barplot.png",width=6,height=4)

################ Bar plot ################
plotting_df <- gather(df,type,percentage,
                      c("deletion","insertion","substitution")) %>%
  group_by(platform,type) %>%
  summarise(sd=sd(percentage),percentage=mean(percentage)) %>%
  ungroup

order <- group_by(plotting_df,platform) %>%
  summarise(percentage=sum(percentage)) %>%
  ungroup %>%
  arrange(percentage) %>%
  pull(platform)
ggplot(plotting_df,aes(x=factor(platform,levels=order),y=percentage,
                       fill=factor(type,levels=c("substitution","insertion","deletion")))) + 
  geom_bar(stat="identity",position="dodge") +
  geom_errorbar(aes(ymin=percentage,ymax=percentage+sd),width=.2,
                position=position_dodge(.9),color="gray") +
  scale_y_continuous(breaks=seq(0,0.6,0.1)) +
  coord_flip() +
  scale_fill_manual(values=c("deletion"="#F8766D",
                             "insertion"="#00BFC4",
                             "substitution"="#7CAE00")) +
  theme(
    # panel.background = element_blank(),
    # panel.grid.minor = element_line(color="black"),
    legend.position="bottom",
    axis.line = element_line(color="black"),
    axis.text = element_text(size=12),
    axis.line.y = element_blank(),
    # axis.text.x = element_text(angle=315),
    axis.title = element_blank(),
    legend.text = element_text(size=20),
    legend.title = element_blank())
ggsave("./figures/fig2_mutation_percentage_barplot.png",width=6,height=4)
