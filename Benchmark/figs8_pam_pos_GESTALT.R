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

df <- data.frame()
for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
    tmp <- GESTALT[[sample_name]] %>%
      filter(type!="unmutated") %>%
      filter(start>2) %>%
      group_by(type,start,length,mutate_to) %>%
      summarise(count=sum(count)) %>%
      ungroup %>%
      mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
             ending_sgRNA=findInterval(end,target_end)+1) %>%
      mutate(target_type=case_when(starting_sgRNA==ending_sgRNA ~ "intra-target",
                                   TRUE ~ "inter-target")) %>%
      mutate(starting_cutsite=start+1-target_end[starting_sgRNA],
             ending_cutsite=end+1-target_end[ending_sgRNA],
             frequency=count/sum(count),
             sample_name=sample_name)
    df <- rbind(df,tmp)}

inter_target <- filter(df,type=="deletion",target_type=="inter-target") %>% 
  group_by(starting_sgRNA,ending_sgRNA) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  mutate(frequency=count/sum(count))

node_size <- gather(inter_target,tmp,target,starting_sgRNA,ending_sgRNA) %>%
  group_by(target) %>%
  summarise(size=sum(frequency)) %>%
  ungroup
col_fun = colorRamp2(c(0, 0.05), c("white", "firebrick2"))
png("./figures/figs8_arc_GESTALT.png",width=600,height=200)
arcplot(as.matrix(select(inter_target,starting_sgRNA,ending_sgRNA)), ordering=1:10,
        lwd.arcs=100*inter_target$frequency, col.arcs=col_fun(inter_target$frequency),
        alpha=0.2,show.nodes=FALSE,show.labels=FALSE,
        col.nodes="gray80",bg.nodes="gray90")
dev.off()


inter <- filter(df,type=="deletion",target_type=="inter-target") %>%
  # filter(length>20) %>%
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
l <- length(unique(inter$starting_cutsite))
mat <- matrix(0,nrow=l,ncol=l)
mat[cbind(inter$starting_cutsite+l,inter$ending_cutsite+l)] <- inter$avg_frequency
split <- factor(c(rep("spacer",20),rep("PAM",3),rep("link",4)),levels=c("spacer","PAM","link"))
png("./figures/figs8_heatmap_inter_GESTALT.png",height=300,width=300)
Heatmap(mat, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE,column_split = split,row_split = split,
        row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),border=TRUE)
dev.off()

png("./figures/fig4_heatmap_legend.png")
l <- Legend(col_fun=col_fun, direction = "horizontal")
draw(l)
dev.off()


intra <- filter(df,type=="deletion",target_type=="intra-target") %>%
  group_by(sample_name) %>%
  mutate(frequency=frequency/sum(frequency)) %>%
  ungroup %>%
  group_by(sample_name,starting_cutsite,ending_cutsite) %>%
  mutate(frequency=sum(frequency)) %>%
  ungroup %>%
  group_by(starting_cutsite,ending_cutsite) %>%
  summarise(avg_frequency=mean(frequency),sd_frequency=sd(frequency)) %>%
  ungroup
mat <- matrix(0,nrow=l,ncol=l)
mat[cbind(intra$starting_cutsite+l,intra$ending_cutsite+l)] <- intra$avg_frequency
png("./figures/figs8_heatmap_intra_GESTALT.png",height=300,width=300)
Heatmap(mat, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE,column_split = split,row_split = split,
        row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),border=TRUE)
dev.off()


inter_start <- filter(df,type=="deletion",target_type=="inter-target") %>%
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
ggplot(inter_start) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(0,0.3)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("./figures/figs8_pam_inter_start_GESTALT.png",width=6,height=1.5)    

inter_end <- filter(df,type=="deletion",target_type=="inter-target") %>%
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
ggplot(inter_end) + 
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
ggsave("./figures/figs8_pam_inter_end_GESTALT.png",width=6,height=1.5)    

intra_start <- filter(df,type=="deletion",target_type=="intra-target") %>%
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

ggplot(intra_start) + 
  geom_bar(aes(x=cutsite,y=avg_frequency),stat="identity") +
  geom_errorbar(aes(x=cutsite,ymin=avg_frequency,ymax=avg_frequency+sd_frequency),width=.3) +
  coord_cartesian(ylim=c(0,0.3)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("./figures/figs8_pam_intra_start_GESTALT.png",width=6,height=1.5)

intra_end <- filter(df,type=="deletion",target_type=="intra-target") %>%
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
ggplot(intra_end) + 
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
ggsave("./figures/figs8_pam_intra_end_GESTALT.png",width=6,height=1.5)  


substitution <- filter(df,type=="substitution") %>%
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
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("./figures/figs8_pam_substitution_GESTALT.png",width=6,height=1.5)  


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
  coord_cartesian(ylim=c(0,0.6)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title = element_blank())
ggsave("./figures/figs8_pam_insertion_GESTALT.png",width=6,height=1.5)  
