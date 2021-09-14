library(ggplot2)
library(dplyr)
library(stringr)
library(TraceQC)
library(tidyr)
library(purrr)

mutated_percentage <- data.frame()
mutation_per_sequence <- data.frame()
mutated_pam_percentage <- data.frame()
ref_length <- 1000
################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_pam_end <- filter(ref$regions,region=="pam") %>%
   pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
   pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
   pull(end)
pam_pos2 <- ref_pam_end-ref_target_start+1
pam_pos1 <- pam_pos2-1

pam_percentage <- function(df) {
   if ("unmutated" %in% df$type) {return(0)} else {
   res <- rep(0,ref_length)
   start <- df$start
   end <- df$start + df$length - 1
   for (i in 1:length(start)) {
      res[start[i]:end[i]] <- 1} 
   return(sum(res[pam_pos1]|res[pam_pos2])/10)}}

for (sample_name in names(Carlin)[2:length(Carlin)]) {
   if (str_detect(sample_name,"0h")|str_detect(sample_name,"G0")) {
      unmutated_count <- Carlin[[sample_name]] %>%
         filter(type=="unmutated") %>%
         pull(count)
      tmp <- Carlin[[sample_name]] %>%
         filter(type!="unmutated") %>%
         group_by(target_seq) %>%
         summarise(mutation_per_sequence=n(),count=mean(count)) %>%
         ungroup
      mutated_percentage <- rbind(mutated_percentage,data.frame(
         mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
                                time="t0",platform="Carlin"))
      mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
         mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
         time="t0",platform="Carlin"))
      
      tmp <- filter(Carlin[[sample_name]]) %>%
         group_by(target_seq) %>%
         nest() %>%
         mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map(data,pam_percentage))) %>%
         ungroup
      mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
         mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
         time="t0",platform="Carlin"))
      }}

for (sample_name in names(Carlin)[2:length(Carlin)]) {
   if ((str_detect(sample_name,"96h")&(!str_detect(sample_name,"NegDox"))|str_detect(sample_name,"G3"))) {
      unmutated_count <- Carlin[[sample_name]] %>%
         filter(type=="unmutated") %>%
         pull(count)
      tmp <- Carlin[[sample_name]] %>%
         filter(type!="unmutated") %>%
         group_by(target_seq) %>%
         summarise(mutation_per_sequence=n(),count=mean(count)) %>%
         ungroup
      mutated_percentage <- rbind(mutated_percentage,data.frame(
         mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
         time="t1",platform="Carlin"))
      mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
         mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
         time="t1",platform="Carlin"))
      tmp <- filter(Carlin[[sample_name]]) %>%
         group_by(target_seq) %>%
         nest() %>%
         mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map(data,pam_percentage))) %>%
         ungroup
      mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
         mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
         time="t1",platform="Carlin"))
      }}

################ hgRNA_invitro ################ 
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
all_pam_pos <- list()
# 103-105
for (ref_file in list.files("../hgRNA-invitro/data/000_ref/")) {
   ref <- parse_ref_file(sprintf("../hgRNA-invitro/data/000_ref/%s",ref_file))
   target_start <- ref$regions %>% filter(.data$region == "target") %>% pull(.data$start)
   target_end <- ref$regions %>% filter(.data$region == "target") %>% pull(.data$end)
   all_pam_pos[[file_path_sans_ext(ref_file)]] <- c(target_end-135+104-target_start+1,target_end-135+104-target_start+2)}
pam_percentage <- function(df,ref) {
   if ("unmutated" %in% df$type) {return(0)} else {
      res <- rep(0,ref_length)
      start <- df$start
      end <- df$start + df$length - 1
      for (i in 1:length(start)) {
         res[start[i]:end[i]] <- 1} 
      return(sum(res[all_pam_pos[[ref]][1]]|res[all_pam_pos[[ref]][2]]))}}

for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
   if (str_detect(sample_name,"0d")&str_detect(sample_name,"A")) {
      unmutated_count <- hgRNA_invitro[[sample_name]] %>%
         filter(type=="unmutated") %>%
         pull(count)
      tmp <- hgRNA_invitro[[sample_name]] %>%
         filter(type!="unmutated") %>%
         group_by(target_seq) %>%
         summarise(mutation_per_sequence=n(),count=mean(count)) %>%
         ungroup
      mutated_percentage <- rbind(mutated_percentage,data.frame(
         mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
         time="t0",platform="hgRNA-invitro"))
      mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
         mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
         time="t0",platform="hgRNA-invitro"))
      ref <- strsplit(sample_name,split="-")[[1]][1]
      tmp <- hgRNA_invitro[[sample_name]] %>%
         group_by(target_seq) %>%
         nest() %>%
         mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map2(data,ref,pam_percentage))) %>%
         ungroup
      mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
         mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
         time="t0",platform="hgRNA-invitro"))
      }}

for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
   if (str_detect(sample_name,"14d")) {
      unmutated_count <- hgRNA_invitro[[sample_name]] %>%
         filter(type=="unmutated") %>%
         pull(count)
      tmp <- hgRNA_invitro[[sample_name]] %>%
         filter(type!="unmutated") %>%
         group_by(target_seq) %>%
         summarise(mutation_per_sequence=n(),count=mean(count)) %>%
         ungroup
      mutated_percentage <- rbind(mutated_percentage,data.frame(
         mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
         time="t1",platform="hgRNA-invitro"))
      mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
         mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
         time="t1",platform="hgRNA-invitro"))
      tmp <- hgRNA_invitro[[sample_name]] %>%
         group_by(target_seq) %>%
         nest() %>%
         mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map2(data,ref,pam_percentage))) %>%
         ungroup
      mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
         mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
         time="t1",platform="hgRNA-invitro"))}}

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
ref_pam_end <- filter(ref$regions,region=="Pam") %>%
   pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
   pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
   pull(end)
pam_pos2 <- ref_pam_end-ref_target_start+1
pam_pos1 <- pam_pos2-1

pam_percentage <- function(df) {
   if ("unmutated" %in% df$type) {return(0)} else {
      res <- rep(0,ref_length)
      start <- df$start
      end <- df$start + df$length - 1
      for (i in 1:length(start)) {
         res[start[i]:end[i]] <- 1} 
      return(sum(res[pam_pos1]|res[pam_pos2],na.rm=TRUE)/5)}}

for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
      unmutated_count <- GESTALT[[sample_name]] %>%
         filter(type=="unmutated") %>%
         pull(count)
      if (length(unmutated_count)==0) {unmutated_count<-0}
      
      tmp <- GESTALT[[sample_name]] %>%
         filter(type!="unmutated") %>%
         group_by(target_seq) %>%
         summarise(mutation_per_sequence=n(),count=mean(count)) %>%
         ungroup
      mutated_percentage <- rbind(mutated_percentage,data.frame(
         mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
         time="t1",platform="GESTALT"))
      mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
         mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
         time="t1",platform="GESTALT"))
      
      tmp <- filter(GESTALT[[sample_name]]) %>%
         group_by(target_seq) %>%
         nest() %>%
         mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map(data,pam_percentage))) %>%
         ungroup
      mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
         mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
         time="t1",platform="GESTALT"))}

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
all_pam_pos <- list()
for (ref_file in list.files("../hgRNA-invivo/data/000_ref/",pattern="L[0-9][0-9]_")) {
   ref <- parse_ref_file(sprintf("../hgRNA-invivo/data/000_ref/%s",ref_file))
   identifier <- file_path_sans_ext(ref_file)
   identifier <- strsplit(identifier,split="_")[[1]][2]
   target_start <- ref$regions %>% filter(.data$region == "target") %>% pull(.data$start)
   spacer_end <- ref$regions %>% filter(.data$region == "spacer") %>% pull(.data$end)
   all_pam_pos[[identifier]] <- c(spacer_end-target_start+3,spacer_end-target_start+4)}

pam_percentage <- function(df) {
   if ("unmutated" %in% df$type) {return(0)} else {
      res <- rep(0,ref_length)
      start <- df$start
      end <- df$start + df$length - 1
      for (i in 1:length(start)) {
         res[start[i]:end[i]] <- 1}
      return(sum(res[all_pam_pos[[df$identifier[1]]][1]]|res[all_pam_pos[[df$identifier[1]]][2]]))}}

info_table <- hgRNA_invivo$info_table %>%
   filter(dev_stage=="Adult")
for (sample_name in info_table$Run) {
   if (sample_name%in%names(hgRNA_invivo)) {
      if (nrow(hgRNA_invivo[[sample_name]])>100) {
   unmutated_count <- hgRNA_invivo[[sample_name]] %>%
      filter(type=="unmutated") %>%
      pull(count) %>%
      sum()
   if (length(unmutated_count)==0) {unmutated_count<-0}
   
   tmp <- hgRNA_invivo[[sample_name]] %>%
      filter(type!="unmutated") %>%
      group_by(target_seq) %>%
      summarise(mutation_per_sequence=n(),count=mean(count)) %>%
      ungroup
   mutated_percentage <- rbind(mutated_percentage,data.frame(
      mutated_percentage=1-unmutated_count/(unmutated_count+sum(tmp$count)),
      time="t1",platform="hgRNA-invivo"))
   mutation_per_sequence <- rbind(mutation_per_sequence,data.frame(
      mutation_per_sequence=sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count),
      time="t1",platform="hgRNA-invivo"))
   tmp <- hgRNA_invivo[[sample_name]] %>%
      group_by(target_seq) %>%
      nest() %>%
      mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map(data,pam_percentage))) %>%
      ungroup
   mutated_pam_percentage <- rbind(mutated_pam_percentage,data.frame(
      mutation_per_sequence=sum(tmp$mutated_pam*tmp$count)/sum(tmp$count),
      time="t1",platform="hgRNA_invivo"))}}}

plotting_df <- group_by(mutated_percentage,platform,time) %>%
   summarise(percentage=mean(mutated_percentage),sd=sd(mutated_percentage)) %>%
   ungroup %>%
   arrange(percentage) %>%
   mutate(inducible=case_when(platform%in%c("Carlin","hgRNA-invitro") ~ "inducible",
                              TRUE ~ "non-inducible")) %>%
   mutate(x=c(0,1,3,4,5,6))
   
ggplot(plotting_df,aes(x=x,y=percentage,fill=inducible)) +
   geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=percentage,ymax=percentage+sd),width=.2) +
   scale_fill_brewer(palette="Paired") +
   scale_x_continuous(breaks=c(0.5,1.5,3.5,4.5,5.5,6.5),labels=plotting_df$platform) +
   geom_vline(xintercept=2,linetype="dashed") +
   theme(
      panel.background = element_blank(),
      axis.line = element_line(color="black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=15),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none")
ggsave("./figures/fig3_mutated_percentage.png",width=3,height=3)

plotting_df <- group_by(mutation_per_sequence,platform,time) %>%
   summarise(avg=mean(mutation_per_sequence),sd=sd(mutation_per_sequence)) %>%
   ungroup %>%
   arrange(avg) %>%
   mutate(inducible=case_when(platform%in%c("Carlin","hgRNA-invitro") ~ "inducible",
                              TRUE ~ "non-inducible")) %>%
   mutate(x=c(0,3,1,4,5,6))
ggplot(plotting_df,aes(x=x,y=avg,fill=inducible)) +
   geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=avg,ymax=avg+sd),width=.2) +
   scale_fill_brewer(palette="Paired") +
   scale_x_continuous(breaks=c(0.5,1.5,3.5,4.5,5.5,6.5),labels=plotting_df$platform) +
   geom_vline(xintercept=2,linetype="dashed") +
   theme(
      panel.background = element_blank(),
      axis.line = element_line(color="black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=15),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none")
ggsave("./figures/fig3_mutation_per_sequence.png",width=3,height=3)

plotting_df <- group_by(mutated_pam_percentage,platform,time) %>%
   summarise(avg=mean(mutation_per_sequence),sd=sd(mutation_per_sequence)) %>%
   ungroup %>%
   arrange(avg) %>%
   mutate(inducible=case_when(platform%in%c("Carlin","hgRNA-invitro") ~ "inducible",
                              TRUE ~ "non-inducible")) %>%
   mutate(x=c(0,1,3,4,5,6))
ggplot(plotting_df,aes(x=x,y=avg,fill=inducible)) +
   geom_bar(stat="identity") +
   geom_errorbar(aes(ymin=avg,ymax=avg+sd),width=.2) +
   scale_fill_brewer(palette="Paired") +
   scale_x_continuous(breaks=c(0.5,1.5,3.5,4.5,5.5,6.5),labels=plotting_df$platform) +
   geom_vline(xintercept=2,linetype="dashed") +
   theme(
      panel.background = element_blank(),
      axis.line = element_line(color="black"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=15),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=20),
      legend.position = "right")
ggsave("./figures/fig3_mutated_pam.png",width=5.5,height=3)
