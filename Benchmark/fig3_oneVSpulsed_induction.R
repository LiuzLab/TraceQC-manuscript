library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(purrr)

ref_length <- 1000
############## Carlin ##############
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_pam_end <- filter(ref$regions,region=="pam") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
ref_length <- ref_target_end - ref_target_start + 1
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

data_summary <- data.frame()
for (sample_name in names(Carlin)[2:length(Carlin)]) {
  sample <- Carlin[[sample_name]] %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(cigar) %>%
    summarise(count=sum(count))
  counts <- sample$count
  shannon_entropy <- Entropy(counts)
  gini_index <- Gini(counts)
  diversity <- length(counts)
  
  unmutated_count <- Carlin[[sample_name]] %>%
    filter(type=="unmutated") %>%
    pull(count)
  tmp <- Carlin[[sample_name]] %>%
    filter(type!="unmutated") %>%
    group_by(target_seq) %>%
    summarise(mutation_per_sequence=n(),count=mean(count)) %>%
    ungroup
  mutated_percentage <- 1-unmutated_count/(unmutated_count+sum(tmp$count))
  mutation_per_sequence <- sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count)
  
  tmp <- Carlin[[sample_name]] %>%
    group_by(target_seq) %>%
    nest() %>%
    mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map(data,pam_percentage))) %>%
    ungroup
  mutated_pam_percentage <- sum(tmp$mutated_pam*tmp$count)/sum(tmp$count)
  
  tmp <- filter(Carlin[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0),
           substitution=replace_na(substitution,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  
  data_summary <- rbind(data_summary,
                        data.frame(sample_name=sample_name,
                                   shannon_entropy=shannon_entropy,
                                   gini_index=gini_index,
                                   diversity=diversity,
                                   avg_deletion_length=avg_deletion_length,
                                   mutated_percentage=mutated_percentage,
                                   mutation_per_sequence=mutation_per_sequence,
                                   mutated_pam_percentage=mutated_pam_percentage))}

Chronic_induction <- filter(data_summary,str_detect(sample_name,"ChronicInduction"))
tmp <- strsplit(Chronic_induction$sample_name,"_")
x_labels <- c(0,3,6,12,18,24)
names(x_labels)=c("0h","12h","24h","48h","72h","96h")
group <- lapply(tmp,function(x) {paste(x[1:(length(x)-1)],collapse="_")}) %>%
  unlist()
time <- lapply(tmp,function(x) {x[length(x)]}) %>%
  unlist()
Chronic_induction$group <- group
Chronic_induction$time <- x_labels[time]

ChronicInduction_PosDox_0h <- filter(Chronic_induction,sample_name=="ChronicInduction_PosDox_0h")
Chronic_induction <- filter(Chronic_induction,group!="ChronicInduction_PosDox") %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_High")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_Low")) %>%
  rbind(mutate(ChronicInduction_PosDox_0h,group="ChronicInduction_PosDox_Medium"))

df <- data_summary %>%
  filter(str_detect(sample_name,"PulsedInduction")) %>%
  separate(sample_name,c("sample_name","time"),"_")
Pulsed_induction <- data.frame()
g0 <- filter(df,time=="G0") %>%
  mutate(time=0)
Pulsed_induction <- rbind(Pulsed_induction,mutate(g0,group="A1")) %>%
  rbind(mutate(g0,group="A2")) %>%
  rbind(mutate(g0,group="B1")) %>%
  rbind(mutate(g0,group="B2"))

g1 <- filter(df,time=="G1") %>%
  mutate(time=8)
Pulsed_induction <- rbind(Pulsed_induction,mutate(g1,group="A1")) %>%
  rbind(mutate(g1,group="A2")) %>%
  rbind(mutate(g1,group="B1")) %>%
  rbind(mutate(g1,group="B2"))

g2a <- filter(df,time=="G2A") %>%
  mutate(time=16)
Pulsed_induction <- rbind(Pulsed_induction,mutate(g2a,group="A1")) %>%
  rbind(mutate(g2a,group="A2"))

g2b <- filter(df,time=="G2B") %>%
  mutate(time=16)
Pulsed_induction <- rbind(Pulsed_induction,mutate(g2b,group="B1")) %>%
  rbind(mutate(g2b,group="B2"))

g3 <- filter(df,str_detect(time,"3")) %>%
  mutate(group=substr(time,start=nchar(time)-1,stop=nchar(time)),
         time=24)
Pulsed_induction <- rbind(Pulsed_induction,g3)
tmp <- Pulsed_induction$time
Pulsed_induction$time <- NULL
Pulsed_induction$time <- tmp
df <- rbind(Chronic_induction,Pulsed_induction) %>%
  filter(!(group %in% c("ChronicInduction_NegDox_Trial1","ChronicInduction_NegDox_Trial2"))) %>%
  mutate(Induction=str_detect(sample_name,"ChronicInduction"))

ggplot(data=df,aes(x=time,y=mutated_percentage,
                   group=group,color=Induction)) +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values=c("TRUE"="#E69F00","FALSE"="#56B4E9")) +
  # annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  scale_x_continuous(breaks=c("t0"=0,"t1"=8,"t2"=16,"t3"=24),
                     sec.axis=dup_axis(breaks=
                                         c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.top = element_text(size=15,color="#E69F00"),
    axis.text.x.bottom = element_text(size=15,color="#56B4E9"),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig3_Carlin_oneVSpulsed_mutated_percent.png",height=3,width=3)

ggplot(data=df,aes(x=time,y=mutated_pam_percentage,
                   group=group,color=Induction)) +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values=c("TRUE"="#E69F00","FALSE"="#56B4E9")) +
  # annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  scale_x_continuous(breaks=c("t0"=0,"t1"=8,"t2"=16,"t3"=24),
                     sec.axis=dup_axis(breaks=
                                         c("0"=0,"12"=3,"24"=6,"48"=12,"72"=18,"96h"=24))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.top = element_text(size=15,color="#E69F00"),
    axis.text.x.bottom = element_text(size=15,color="#56B4E9"),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig3_Carlin_oneVSpulsed_mutated_PAM.png",height=3,width=3)

############## hgRNA-invitro ##############
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
all_pam_pos <- list()
ref_length <- 1000

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

data_summary <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  tmp <- strsplit(sample_name,"-") %>% unlist()
  barcode <- tmp[1]
  day <- tmp[2]
  if (str_detect(sample_name,"\\(A'\\)")) {
    barcode <- "A'21"
    day <- substr(day,start=1,stop=nchar(day)-4)}
  sample <- hgRNA_invitro[[sample_name]] %>%
    mutate(cigar=paste(type,start,length,mutate_to,sep="_")) %>%
    group_by(cigar) %>%
    summarise(count=sum(count))
  counts <- sample$count
  shannon_entropy <- Entropy(counts)
  gini_index <- Gini(counts)
  diversity <- length(counts)
  
  unmutated_count <- hgRNA_invitro[[sample_name]] %>%
    filter(type=="unmutated") %>%
    pull(count)
  tmp <- hgRNA_invitro[[sample_name]] %>%
    filter(type!="unmutated") %>%
    group_by(target_seq) %>%
    summarise(mutation_per_sequence=n(),count=mean(count)) %>%
    ungroup
  mutated_percentage <- 1-unmutated_count/(unmutated_count+sum(tmp$count))
  if (length(mutated_percentage)==0) {mutated_percentage<-0}
  mutation_per_sequence <- sum(tmp$mutation_per_sequence*tmp$count)/sum(tmp$count)
  ref <- strsplit(sample_name,split="-")[[1]][1]
  tmp <- filter(hgRNA_invitro[[sample_name]]) %>%
    group_by(target_seq) %>%
    nest() %>%
    mutate(count=unlist(map(data,function(x) {mean(x$count)})),mutated_pam=unlist(map2(data,ref,pam_percentage))) %>%
    ungroup
  mutated_pam_percentage <- sum(tmp$mutated_pam*tmp$count)/sum(tmp$count)
  
  tmp <- filter(hgRNA_invitro[[sample_name]],type!="unmutated") %>%
    group_by(target_seq,type) %>%
    summarise(length=sum(length),count=mean(count)) %>%
    ungroup %>%
    spread(type,length) %>%
    mutate(deletion=replace_na(deletion,0),
           insertion=replace_na(insertion,0))
  avg_deletion_length <- sum(tmp$deletion*tmp$count)/sum(tmp$count)
  
  data_summary <- rbind(data_summary,
                        data.frame(sample_name=sample_name,
                                   barcode=barcode,time=day,
                                   shannon_entropy=shannon_entropy,
                                   gini_index=gini_index,
                                   diversity=diversity,
                                   avg_deletion_length=avg_deletion_length,
                                   mutated_percentage=mutated_percentage,
                                   mutation_per_sequence=mutation_per_sequence,
                                   mutated_pam_percentage=mutated_pam_percentage))}

pulsed_induction <- filter(data_summary,barcode=="A21",str_detect(time,"pop"))
pop3 <- filter(pulsed_induction,time=="pop3")
pop3 <- rbind(mutate(pop3,time=0,group="pop321"),
              mutate(pop3,time=0,group="pop322"),
              mutate(pop3,time=0,group="pop323"))
pop32 <- filter(pulsed_induction,time=="pop32")
pop32 <- rbind(mutate(pop32,time=7,group="pop321"),
               mutate(pop32,time=7,group="pop322"),
               mutate(pop32,time=7,group="pop323"))
df <- filter(pulsed_induction,str_detect(time,"pop32[:digit:]")) %>%
  mutate(group=time,time=14) %>%
  rbind(pop3,pop32)

pop4 <- filter(pulsed_induction,time=="pop4")
pop4 <- rbind(mutate(pop3,time=0,group="pop411"),
              mutate(pop3,time=0,group="pop412"),
              mutate(pop3,time=0,group="pop413"))
pop41 <- filter(pulsed_induction,time=="pop41")
pop41 <- rbind(mutate(pop41,time=7,group="pop411"),
               mutate(pop41,time=7,group="pop412"),
               mutate(pop41,time=7,group="pop413"))
df <- filter(pulsed_induction,str_detect(time,"pop41[:digit:]")) %>%
  mutate(group=time,time=14) %>%
  rbind(pop4,pop41,df)

x_labels <- c(0,2,14)
names(x_labels)=c("0d","2d","14d")
df2 <- filter(data_summary,barcode=="A\'21") %>%
  mutate(time=x_labels[time],group="one_time")
df <- rbind(df,df2) %>%
  mutate(Induction=case_when(str_detect(group,"pop") ~ "pulsed",
                             TRUE ~ "one-time"))
  
ggplot(data=df,aes(x=time,y=mutated_percentage,group=group,color=Induction)) +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values=c("pulsed"="#56B4E9","one-time"="#E69F00")) +
  scale_x_continuous(breaks=c("t0"=0,"t1"=7,"t2"=14),
                     sec.axis=dup_axis(breaks=
                                         c("0d"=0,"2d"=2,"14d"=14))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.top = element_text(size=15,color="#E69F00"),
    axis.text.x.bottom = element_text(size=15,color="#56B4E9"),
    axis.title = element_blank(),
    legend.position = "none")
ggsave("./figures/fig3_hgRNA-invitro_oneVSpulsed_mutated_percent.png",height=3,width=3)

ggplot(data=df,aes(x=time,y=mutated_pam_percentage,group=group,color=Induction)) +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values=c("pulsed"="#56B4E9","one-time"="#E69F00")) +
  scale_x_continuous(breaks=c("t0"=0,"t1"=7,"t2"=14),
                     sec.axis=dup_axis(breaks=
                                         c("0d"=0,"2d"=2,"14d"=14))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.top = element_text(size=15,color="#E69F00"),
    axis.text.x.bottom = element_text(size=15,color="#56B4E9"),
    axis.title = element_blank(),
    legend.position = "bottom")
ggsave("./figures/fig3_hgRNA-invitro_oneVSpulsed_mutated_pam_percent.png",height=6,width=3)

ggplot(data=df,aes(x=time,y=mutated_pam_percentage,group=group,color=Induction)) +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values=c("pulsed"="#56B4E9","one-time"="#E69F00")) +
  scale_x_continuous(breaks=c("t0"=0,"t1"=7,"t2"=14),
                     sec.axis=dup_axis(breaks=
                                         c("0d"=0,"2d"=2,"14d"=14))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x.top = element_text(size=15,color="#E69F00"),
    axis.text.x.bottom = element_text(size=15,color="#56B4E9"),
    axis.title = element_blank(),
    legend.position = "bottom")
ggsave("./figures/fig3_hgRNA-invitro_oneVSpulsed_mutated_pam_percent.png",height=6,width=3)
