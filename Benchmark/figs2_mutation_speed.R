library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)


############## GESTALT v1 ##############
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_cell_culture.rds")
data_summary <- data.frame()
ref <- parse_ref_file("../GESTALT/data/000_ref_v1.txt")
ref_link_end <- filter(ref$regions,region=="Link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2

target_mutation <- function(data) {
  is_mutated <- rep(0,10)
  for (i in 1:nrow(data)) {
    is_mutated[data$starting_target[i]:data$ending_target[i]] <- 1}
  data.frame(target=1:10,is_mutated=is_mutated)}

for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  df <- GESTALT[[sample_name]] %>%
    filter(type!="unmutated") %>%
    filter(type!="substitution") %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    group_by(target_seq) %>%
    nest() %>%
    mutate(target_array_mutation=map(data,target_mutation)) %>%
    mutate(count=map_dbl(data,function(data) {mean(data$count)})) %>%
    ungroup %>%
    unnest(target_array_mutation) %>%
    group_by(target) %>%
    summarise(mutated_count=sum(count*is_mutated),total_count=sum(count)) %>%
    ungroup
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}

df <- filter(data_summary,sample_name!="SRR3561162") %>%
  mutate(frequency=mutated_count/total_count) %>%
  group_by(target) %>%
  summarise(`mutated percentage`=mean(frequency),sd=sd(frequency)) %>%
  ungroup

ggplot(df,aes(x=as.factor(target),y=`mutated percentage`)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=`mutated percentage`,ymax=`mutated percentage`+sd),width=.2) +
  xlab("targets") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15))
ggsave("./figures/figs2_GESTALT_v6_mutated_targets.png",width=5,height=3)

############## GESTALT v6 ##############
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

target_mutation <- function(data) {
  is_mutated <- rep(0,10)
  for (i in 1:nrow(data)) {
    is_mutated[data$starting_target[i]:data$ending_target[i]] <- 1}
  data.frame(target=1:10,is_mutated=is_mutated)}

for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  df <- GESTALT[[sample_name]] %>%
    filter(type!="unmutated") %>%
    # filter(type!="substitution") %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    group_by(target_seq) %>%
    nest() %>%
    mutate(target_array_mutation=map(data,target_mutation)) %>%
    mutate(count=map_dbl(data,function(data) {mean(data$count)})) %>%
    ungroup %>%
    unnest(target_array_mutation) %>%
    group_by(target) %>%
    summarise(mutated_count=sum(count*is_mutated),total_count=sum(count)) %>%
    ungroup
    
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}

data_summary <- mutate(data_summary,frequency=mutated_count/total_count) %>%
  group_by(target) %>%
  summarise(`mutated percentage`=mean(frequency),sd=sd(frequency)) %>%
  ungroup

ggplot(data_summary,aes(x=as.factor(target),y=`mutated percentage`)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=`mutated percentage`,ymax=`mutated percentage`+sd),width=.2) +
  xlab("targets") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15))
ggsave("./figures/figs2_GESTALT_v6_mutated_targets.png",width=5,height=3)

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

target_mutation <- function(data) {
  is_mutated <- rep(0,10)
  for (i in 1:nrow(data)) {
    is_mutated[data$starting_target[i]:data$ending_target[i]] <- 1}
  data.frame(target=1:10,is_mutated=is_mutated)}

for (sample_name in names(Carlin)[2:length(Carlin)]) {
  df <- Carlin[[sample_name]] %>%
    filter(type!="unmutated") %>%
    # filter(type!="substitution") %>%
    mutate(end=start+length-1,starting_target=findInterval(start,target_end)+1,
           ending_target=findInterval(end,target_end)+1) %>%
    mutate(ending_target=case_when(ending_target>10 ~ 10,
                                   TRUE ~ ending_target)) %>%
    group_by(target_seq) %>%
    nest() %>%
    mutate(target_array_mutation=map(data,target_mutation)) %>%
    mutate(count=map_dbl(data,function(data) {mean(data$count)})) %>%
    ungroup %>%
    unnest(target_array_mutation) %>%
    group_by(target) %>%
    summarise(mutated_count=sum(count*is_mutated),total_count=sum(count)) %>%
    ungroup
  
  data_summary <- rbind(data_summary,
                        mutate(df,sample_name=sample_name))}

saturated_df <- filter(data_summary,str_detect(sample_name,"G3")|
                      str_detect(sample_name,"96h")&str_detect(sample_name,"PosDox")) %>%
  mutate(frequency=mutated_count/total_count) %>%
  group_by(target) %>%
  summarise(`mutated percentage`=mean(frequency),sd=sd(frequency)) %>%
  ungroup

ggplot(saturated_df,aes(x=as.factor(target),y=`mutated percentage`)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=`mutated percentage`,ymax=`mutated percentage`+sd),width=.2) +
  xlab("targets") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15))
ggsave("./figures/figs2_Carlin_saturated_mutated_targets.png",width=5,height=3)


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
  rbind(mutate(ChronicInduction_PosDox_0h,treatment="ChronicInduction_PosDox_Medium")) %>%
  mutate(`mutated percentage`=mutated_count/total_count)

x_labels <- c(0,1,2,4,6,8)
names(x_labels)=c("0h","12h","24h","48h","72h","96h")
df <- mutate(Chronic_induction,group=target,time_int=x_labels[time]) %>%
  mutate(target=as.factor(target))
ggplot(data=filter(df,treatment=="ChronicInduction_PosDox_Low"),
       aes(x=time_int,y=`mutated percentage`,color=target,fill=target)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave("./figures/figs2_low_mutation_speed.png",height=5,width=5)

ggplot(data=filter(df,treatment=="ChronicInduction_PosDox_Medium"),
       aes(x=time_int,y=`mutated percentage`,color=target,fill=target)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave("./figures/figs2_medium_mutation_speed.png",height=5,width=5)

ggplot(data=filter(df,treatment=="ChronicInduction_PosDox_High"),
       aes(x=time_int,y=`mutated percentage`,color=target,fill=target)) +
  geom_point(size=5) +
  geom_line() +
  scale_x_continuous(breaks=x_labels) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave("./figures/figs2_high_mutation_speed.png",height=5,width=5)


