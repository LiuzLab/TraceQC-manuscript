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
    filter(type=="deletion") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((starting_target==ending_target) ~ "intra-target deletion",
                          (starting_target!=ending_target) ~ "inter-target deletion"))
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}


data_summary <- group_by(data_summary,type,length) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  group_by(type) %>%
  mutate(frequency=count/sum(count)) %>%
  ungroup %>%
  mutate(is_intra=2*(type=="intra-target deletion")-1) %>%
  mutate(frequency=frequency*is_intra)

ggplot(data_summary,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim=c(-0.3,0.3)) +
  scale_fill_brewer(palette = "Set2") +
  # scale_x_continuous(breaks=seq(0.5,
  #                               max(data_points$length)+0.5),
  #                    labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank())
ggsave("./figures/figs1_del_length_Carlin.png",width=5,height=3)

############## scCarlin ##############
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
data_summary <- data.frame()
ref <- parse_ref_file("../scCarlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2

for (sample_name in names(scCarlin)[2:length(scCarlin)]) {
  df <- scCarlin[[sample_name]] %>%
    filter(type=="deletion") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=n()) %>%
    ungroup %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((starting_target==ending_target) ~ "intra-target deletion",
                          (starting_target!=ending_target) ~ "inter-target deletion"))
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}


data_summary <- group_by(data_summary,type,length) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  group_by(type) %>%
  mutate(frequency=count/sum(count)) %>%
  ungroup %>%
  mutate(is_intra=2*(type=="intra-target deletion")-1) %>%
  mutate(frequency=frequency*is_intra)

ggplot(data_summary,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim=c(-0.3,0.3)) +
  scale_fill_brewer(palette = "Set2") +
  # scale_x_continuous(breaks=seq(0.5,
  #                               max(data_points$length)+0.5),
  #                    labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank())
ggsave("./figures/figs1_del_length_scCarlin.png",width=5,height=3)

############## GESTALT ##############
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
data_summary <- data.frame()
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
ref_link_end <- filter(ref$regions,region=="Link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2
target_end <- target_end[2:11]

for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  df <- GESTALT[[sample_name]] %>%
    filter(type=="deletion") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=sum(count)) %>%
    ungroup %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((starting_target==ending_target) ~ "intra-target deletion",
                          (starting_target!=ending_target) ~ "inter-target deletion"))
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}


data_summary <- group_by(data_summary,type,length) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  group_by(type) %>%
  mutate(frequency=count/sum(count)) %>%
  ungroup %>%
  mutate(is_intra=2*(type=="intra-target deletion")-1) %>%
  mutate(frequency=frequency*is_intra)

ggplot(data_summary,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim=c(-0.3,0.3)) +
  scale_fill_brewer(palette = "Set2") +
  # scale_x_continuous(breaks=seq(0.5,
  #                               max(data_points$length)+0.5),
  #                    labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank())
ggsave("./figures/figs1_del_length_GESTALT.png",width=5,height=3)

############## scGESTALT ##############
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
data_summary <- data.frame()
ref <- parse_ref_file("../scGESTALT/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- ref_link_end-ref_target_start+2
target_end <- target_end[2:10]

for (sample_name in names(scGESTALT)[2:length(scGESTALT)]) {
  df <- scGESTALT[[sample_name]] %>%
    filter(type=="deletion") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=n()) %>%
    ungroup %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((starting_target==ending_target) ~ "intra-target deletion",
                          (starting_target!=ending_target) ~ "inter-target deletion"))
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}


data_summary <- group_by(data_summary,type,length) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  group_by(type) %>%
  mutate(frequency=count/sum(count)) %>%
  ungroup %>%
  mutate(is_intra=2*(type=="intra-target deletion")-1) %>%
  mutate(frequency=frequency*is_intra)

ggplot(data_summary,aes(x=length,y=frequency,fill=type)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim=c(-0.3,0.3)) +
  scale_fill_brewer(palette = "Set2") +
  # scale_x_continuous(breaks=seq(0.5,
  #                               max(data_points$length)+0.5),
  #                    labels = labels[1:(max(data_points$length)+1)]) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=330,vjust=0.6),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_blank())
ggsave("./figures/figs1_del_length_scGESTALT.png",width=5,height=3)
