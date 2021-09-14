library(TraceQC)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(arcdiagram)
library(stringr)
library(arcdiagram)

hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")

all_pam_end_pos <- list()
for (ref_file in list.files("../hgRNA-invitro/data/000_ref/")) {
  ref <- parse_ref_file(sprintf("../hgRNA-invitro/data/000_ref/%s",ref_file))
  target_start <- ref$regions %>% filter(.data$region == "target") %>% pull(.data$start)
  target_end <- ref$regions %>% filter(.data$region == "target") %>% pull(.data$end)
  all_pam_end_pos[[file_path_sans_ext(ref_file)]] <- target_end-target_start-29}
tmp <-substr(ref$refseq,start=target_start,stop=target_end)

df <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- strsplit(sample_name,split="-")[[1]]
  target_end <- all_pam_end_pos[[tmp[1]]]
  if (substr(tmp[1],start=2,stop=3)=="21"&(!str_detect(tmp[2],"pop"))) {
    # print(sample_name)
    # if (tmp[1]=="A21") {
    tmp <- hgRNA_invitro[[sample_name]] %>%
      filter(type!="unmutated") %>%
      filter(start>2) %>%
      group_by(type,start,length,mutate_to) %>%
      summarise(count=sum(count)) %>%
      ungroup %>%
      mutate(end=start+length-1) %>%
      mutate(starting_cutsite=start-target_end,
             ending_cutsite=end-target_end,
             frequency=count/sum(count),
             sample_name=sample_name) %>%
      filter(starting_cutsite<=15,starting_cutsite>=-23) %>%
      filter(ending_cutsite<=15,ending_cutsite>=-23)
    df <- rbind(df,tmp)}}

deletion <- filter(df,type=="deletion") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  group_by(sample_name,starting_cutsite,ending_cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(starting_cutsite,ending_cutsite) %>%
  summarise(avg_frequency=mean(frequency),sd_frequency=sd(frequency)) %>%
  ungroup
col_fun = colorRamp2(c(0, 0.05), c("white", "red"))
l <- length(unique(c(deletion$starting_cutsite,deletion$ending_cutsite)))
mat <- matrix(0,nrow=l,ncol=l)
mat[cbind(deletion$starting_cutsite+24,deletion$ending_cutsite+24)] <- deletion$avg_frequency
split <- factor(c(rep("spacer",21),rep("PAM",3),rep("scaffold",15)),levels=c("spacer","PAM","scaffold"))
png("./figures/figs9_heatmap_del_hgRNA-invitro.png",height=300,width=300)
Heatmap(mat, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE,column_split = split,row_split = split,
        row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),border=TRUE)
dev.off()


del_start <- filter(df,type=="deletion") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  mutate(cutsite=starting_cutsite) %>%
  group_by(sample_name,cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(cutsite) %>%
  summarise(avg_frequency=mean(frequency),sd_frequency=sd(frequency)) %>%
  ungroup
ggplot(del_start) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(0,0.3)) +
  scale_x_continuous(breaks = c(-22,-17,-12,-7,-3,-2,-1,0,5,10,15),
                     labels = c(20,15,10,5,1,"G","G","G",5,10,15)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_blank())
ggsave("./figures/figs9_pam_del_start_hgRNA_invitro.png",width=6,height=1.5)    

del_end <- filter(df,type=="deletion") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  mutate(cutsite=ending_cutsite) %>%
  group_by(sample_name,cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(cutsite) %>%
  summarise(avg_frequency=-mean(frequency),sd_frequency=-sd(frequency)) %>%
  ungroup
ggplot(del_end) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(-0.3,0)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("./figures/figs9_pam_del_end_hgRNA_invitro.png",width=6,height=1.5)    


substitution <- filter(df,type=="mutation") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  mutate(cutsite=starting_cutsite) %>%
  group_by(sample_name,cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(cutsite) %>%
  summarise(avg_frequency=mean(frequency),sd_frequency=sd(frequency)) %>%
  ungroup
ggplot(substitution) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(0,0.3)) +
  scale_x_continuous(breaks = c(-22,-17,-12,-7,-3,-2,-1,0,5,10,15),
                     labels = c(20,15,10,5,1,"G","G","G",5,10,15)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_blank())
ggsave("./figures/figs9_pam_substitution_hgRNA_invitro.png",width=6,height=1.5)  


insertion <- filter(df,type=="insertion") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  mutate(cutsite=starting_cutsite) %>%
  group_by(sample_name,cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(cutsite) %>%
  summarise(avg_frequency=mean(frequency),sd_frequency=sd(frequency)) %>%
  ungroup
ggplot(insertion) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_continuous(breaks = c(-22,-17,-12,-7,-3,-2,-1,0,5,10,15),
                     labels = c(20,15,10,5,1,"G","G","G",5,10,15)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_blank())
ggsave("./figures/figs9_pam_insertion_hgRNA_invitro.png",width=6,height=1.5)  
