res <- data.frame()

############## Carlin ##############
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2

all_samples <- list()
for (sample_name in names(Carlin)[2:length(names(Carlin))]) {
  all_samples[[sample_name]] <- Carlin[[sample_name]] %>%
    filter(type!="unmutated")  %>%
    mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
           ending_sgRNA=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                          (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                          TRUE ~ type)) %>%
    mutate(cigar=paste(type,start,length,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise(log10_count=log10(sum(count)+1)) %>%
    ungroup %>%
    group_by(type) %>%
    arrange(desc(log10_count)) %>%
    mutate(tmp=1:n()) %>%
    ungroup}

samples <- names(all_samples)
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {
    tmp <- strsplit(samples[i],split="_")[[1]]
    s1 <- paste(tmp[1:(length(tmp)-1)],collapse="_")
    tmp <- strsplit(samples[j],split="_")[[1]]
    s2 <- paste(tmp[1:(length(tmp)-1)],collapse="_")
    if (s1 != s2) {
      df1 <- all_samples[[samples[i]]]
      df2 <- all_samples[[samples[j]]]
      df <- full_join(df1,df2,by=c("type","cigar")) %>%
        mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
        mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
        mutate(missing=(log10_count.x==0|log10_count.y==0))
      overlap <- group_by(df,type) %>%
        summarise(value=1-sum(missing)/n()) %>%
        ungroup %>%
        mutate(sample1=samples[i],samples2=samples[j],score="overlap",platform="Carlin")
      pearson <- filter(df,!missing) %>%
        group_by(type) %>%
        summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
        ungroup %>%
        mutate(sample1=samples[i],samples2=samples[j],score="pearson",platform="Carlin")
      res <- rbind(res,overlap,pearson)}}}

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
ref <- parse_ref_file("../scCarlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- ref_link_end-ref_target_start+2

all_samples <- list()
for (sample_name in names(scCarlin)[2:length(names(scCarlin))]) {
  all_samples[[sample_name]] <- scCarlin[[sample_name]] %>%
    filter(type!="unmutated")  %>%
    mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
           ending_sgRNA=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                          (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                          .data$type=="mutation" ~ "substitution",
                          TRUE ~ type)) %>%
    mutate(cigar=paste(type,start,length,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise(log10_count=log10(n()+1)) %>%
    ungroup %>%
    group_by(type) %>%
    arrange(desc(log10_count)) %>%
    mutate(tmp=1:n()) %>%
    ungroup}

samples <- names(all_samples)
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {
    tmp <- strsplit(samples[i],split="_")[[1]]
    s1 <- tmp[1]
    tmp <- strsplit(samples[j],split="_")[[1]]
    s2 <- tmp[1]
    if (s1 != s2) {
      df1 <- all_samples[[samples[i]]]
      df2 <- all_samples[[samples[j]]]
      df <- full_join(df1,df2,by=c("type","cigar")) %>%
        mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
        mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
        mutate(missing=(log10_count.x==0|log10_count.y==0))
      overlap <- group_by(df,type) %>%
        summarise(value=1-sum(missing)/n()) %>%
        ungroup %>%
        mutate(sample1=samples[i],samples2=samples[j],score="overlap",platform="scCarlin")
      pearson <- filter(df,!missing) %>%
        group_by(type) %>%
        summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
        ungroup %>%
        mutate(sample1=samples[i],samples2=samples[j],score="pearson",platform="scCarlin")
      res <- rbind(res,overlap,pearson)}}}

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
ref_link_end <- filter(ref$regions,region=="Link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- c(ref_link_end)-ref_target_start+2
target_end <- target_end[2:11]

all_samples <- list()
for (sample_name in names(GESTALT)[2:length(names(GESTALT))]) {
  all_samples[[sample_name]] <- GESTALT[[sample_name]] %>%
    filter(type!="unmutated")  %>%
    mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
           ending_sgRNA=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                          (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                          .data$type=="mutation" ~ "substitution",
                          TRUE ~ type)) %>%
    mutate(cigar=paste(type,start,length,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise(log10_count=log10(sum(count)+1)) %>%
    ungroup %>%
    group_by(type) %>%
    arrange(desc(log10_count)) %>%
    mutate(tmp=1:n()) %>%
    ungroup}

runs <- names(all_samples)
info <- filter(GESTALT$info_table,Run %in% names(GESTALT))
for (i in 1:(length(runs)-1)) {
  for (j in (i+1):length(runs)) {
    tmp <- filter(info,Run==runs[i]) %>% pull(Experiment)
    s1 <- strsplit(tmp,split="_")[[1]][2]
    tmp <- filter(info,Run==runs[j]) %>% pull(Experiment)
    s2 <- strsplit(tmp,split="_")[[1]][2]
    if (s1 != s2) {
      df1 <- all_samples[[runs[i]]]
      df2 <- all_samples[[runs[j]]]
      df <- full_join(df1,df2,by=c("type","cigar")) %>%
        mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
        mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
        mutate(missing=(log10_count.x==0|log10_count.y==0))
      overlap <- group_by(df,type) %>%
        summarise(value=1-sum(missing)/n()) %>%
        ungroup %>%
        mutate(sample1=runs[i],samples2=runs[j],score="overlap",platform="GESTALT")
      pearson <- filter(df,!missing) %>%
        group_by(type) %>%
        summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
        ungroup %>%
        mutate(sample1=runs[i],samples2=runs[j],score="pearson",platform="GESTALT")
      res <- rbind(res,overlap,pearson)}}}

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
ref <- parse_ref_file("../scGESTALT/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
target_end <- ref_link_end-ref_target_start+2

all_samples <- list()
for (sample_name in names(scGESTALT)[2:length(names(scGESTALT))]) {
  all_samples[[sample_name]] <- scGESTALT[[sample_name]] %>%
    filter(type!="unmutated")  %>%
    mutate(end=start+length-1,starting_sgRNA=findInterval(start,target_end)+1,
           ending_sgRNA=findInterval(end,target_end)+1) %>%
    mutate(type=case_when((.data$type=="deletion"&.data$starting_sgRNA==.data$ending_sgRNA) ~ "intra-target deletion",
                          (.data$type=="deletion"&.data$starting_sgRNA!=.data$ending_sgRNA) ~ "inter-target deletion",
                          .data$type=="mutation" ~ "substitution",
                          TRUE ~ type)) %>%
    mutate(cigar=paste(type,start,length,sep="_")) %>%
    group_by(type,cigar) %>%
    summarise(log10_count=log10(n()+1)) %>%
    ungroup %>%
    group_by(type) %>%
    arrange(desc(log10_count)) %>%
    mutate(tmp=1:n()) %>%
    ungroup}

runs <- names(all_samples)
for (i in 1:(length(runs)-1)) {
  for (j in (i+1):length(runs)) {
    df1 <- all_samples[[runs[i]]]
    df2 <- all_samples[[runs[j]]]
    df <- full_join(df1,df2,by=c("type","cigar")) %>%
      mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
      mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
      mutate(missing=(log10_count.x==0|log10_count.y==0))
    overlap <- group_by(df,type) %>%
      summarise(value=1-sum(missing)/n()) %>%
      ungroup %>%
      mutate(sample1=runs[i],samples2=runs[j],score="overlap",platform="scGESTALT")
    pearson <- filter(df,!missing) %>%
      group_by(type) %>%
      summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
      ungroup %>%
      mutate(sample1=runs[i],samples2=runs[j],score="pearson",platform="scGESTALT")
    res <- rbind(res,overlap,pearson)}}

# ################ hgRNA-invitro ################
# hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
# all_samples <- list()
# for (sample_name in names(hgRNA_invitro)[2:length(names(hgRNA_invitro))]) {
#   if (str_detect(sample_name,"A21")) {
#   all_samples[[sample_name]] <- hgRNA_invitro[[sample_name]] %>%
#     filter(type!="unmutated")  %>%
#     mutate(type=case_when(.data$type=="deletion" ~ "intra-target deletion",
#                           .data$type=="mutation" ~ "substitution",
#                           TRUE ~ type)) %>%
#     mutate(cigar=paste(type,start,length,sep="_")) %>%
#     group_by(type,cigar) %>%
#     summarise(log10_count=log10(sum(count)+1)) %>%
#     ungroup %>%
#     group_by(type) %>%
#     arrange(desc(log10_count)) %>%
#     mutate(tmp=1:n()) %>%
#     ungroup %>%
#     filter(tmp <= 20)}}
# 
# samples <- names(all_samples)
# for (i in 1:(length(samples)-1)) {
#   for (j in (i+1):length(samples)) {
#     if (str_detect(samples[i],"pop")) {
#       s1 <- 1} else if (str_detect(samples[i],"\\(A'\\)")) {
#         s1 <- 2} else {s1 <- 3}
#     if (str_detect(samples[j],"pop")) {
#       s2 <- 1} else if (str_detect(samples[j],"\\(A'\\)")) {
#         s2 <- 2} else {s2 <- 3} 
#     
#     if (s1 != s2) {
#       df1 <- all_samples[[samples[i]]]
#       df2 <- all_samples[[samples[j]]]
#       df <- full_join(df1,df2,by=c("type","cigar")) %>%
#         mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
#         mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
#         mutate(missing=(log10_count.x==0|log10_count.y==0))
#       overlap <- group_by(df,type) %>%
#         summarise(value=1-sum(missing)/n()) %>%
#         ungroup %>%
#         mutate(sample1=samples[i],samples2=samples[j],score="overlap",platform="hgRNA-invitro")
#       pearson <- filter(df,!missing) %>%
#         group_by(type) %>%
#         summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
#         ungroup %>%
#         mutate(sample1=samples[i],samples2=samples[j],score="pearson",platform="hgRNA-invitro")
#       res <- rbind(res,overlap,pearson)}}}
# 
# ################ hgRNA-invivo ################
# hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
# 
# all_samples <- list()
# for (sample_name in names(hgRNA_invivo)[2:length(names(hgRNA_invivo))]) {
#   if (nrow(hgRNA_invivo[[sample_name]])>100) {
#     exp <- filter(hgRNA_invivo$info_table,Run==sample_name) %>% pull(`Sample Name`)
#     if (str_detect(exp,"emb")) {
#       all_samples[[sample_name]] <- hgRNA_invivo[[sample_name]] %>%
#         filter(type!="unmutated")  %>%
#         mutate(type=case_when(.data$type=="deletion" ~ "intra-target deletion",
#                               .data$type=="mutation" ~ "substitution",
#                               TRUE ~ type)) %>%
#         mutate(cigar=paste(type,start,length,sep="_")) %>%
#         group_by(identifier,type,cigar) %>%
#         summarise(log10_count=log10(sum(count)+1)) %>%
#         ungroup}}}
# 
# info <- filter(hgRNA_invivo$info_table,Run %in% names(all_samples))
# samples <- names(all_samples)
# for (i in 1:(length(samples)-1)) {
#   for (j in (i+1):length(samples)) {
#     tmp <- filter(hgRNA_invivo$info_table,Run==samples[i]) %>% pull(`Sample Name`)
#     s1 <- str_extract(tmp,"emb[0-9]")
#     tmp <- filter(hgRNA_invivo$info_table,Run==samples[j]) %>% pull(`Sample Name`)
#     s2 <- str_extract(tmp,"emb[0-9]")
#     
#     if (s1 != s2) {
#       df1 <- all_samples[[samples[i]]]
#       df2 <- all_samples[[samples[j]]]
#       df <- full_join(df1,df2,by=c("identifier","type","cigar")) %>%
#         mutate(log10_count.x=replace_na(log10_count.x,0)) %>%
#         mutate(log10_count.y=replace_na(log10_count.y,0)) %>%
#         mutate(missing=(log10_count.x==0|log10_count.y==0))
#       overlap <- group_by(df,type) %>%
#         summarise(value=1-sum(missing)/n()) %>%
#         ungroup %>%
#         mutate(sample1=samples[i],samples2=samples[j],score="overlap",platform="hgRNA-invivo")
#       pearson <- filter(df,!missing) %>%
#         group_by(type) %>%
#         summarise(value=cor(log10_count.x,log10_count.y,method="pearson")) %>%
#         ungroup %>%
#         mutate(sample1=samples[i],samples2=samples[j],score="pearson",platform="hgRNA-invivo")
#       res <- rbind(res,overlap,pearson)}}}
# 
res2 <- filter(res,!is.na(value),score=="overlap",!str_detect(platform,"hgRNA")) %>%
  filter(type != "substitution")
plotting_df <- group_by(res2,platform,type) %>%
  summarise(avg=mean(value),sd=sd(value)) %>%
  ungroup %>%
  mutate(platform=factor(platform,levels=c("Carlin","scCarlin","GESTALT","scGESTALT")))
p_values <- spread(res2,type,value) %>%
  group_by(platform) %>%
  summarise(p=t.test(.data$`inter-target deletion`,.data$`intra-target deletion`)$p.value) %>%
  ungroup

ggplot(plotting_df,aes(x=platform,y=avg,fill=type)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=avg,ymax=avg+sd),position=position_dodge(width = 0.9),
                width=.2) +
  scale_fill_manual(values=c("inter-target deletion"="gold1","intra-target deletion"="coral4",
                              "insertion"="#00BFC4","substitution"="#7CAE00")) +
  scale_y_continuous(breaks=c(0,0.2,0.4)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15),
    axis.line.x = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
    legend.position = "none")
ggsave("./figures/fig5_overlap.png",width=5,height=3)

res2 <- filter(res,!is.na(value),score=="pearson",!str_detect(platform,"hgRNA"))
plotting_df <- group_by(res2,platform,type) %>%
  summarise(avg=mean(value),sd=sd(value)) %>%
  ungroup %>%
  mutate(platform=factor(platform,levels=c("Carlin","scCarlin","GESTALT","scGESTALT")))
p_values <- spread(res2,type,value) %>%
  group_by(platform) %>%
  summarise(p=t.test(.data$`inter-target deletion`,.data$`intra-target deletion`)$p.value) %>%
  ungroup

ggplot(plotting_df,aes(x=platform,y=avg,fill=type)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=avg,ymax=avg+sd),position=position_dodge(width = 0.9),
                width=.2) +
  # scale_fill_manual(values=c("inter-target deletion"="indianred2",
  #                            "intra-target deletion"="lightskyblue")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15),
    axis.line.x = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
    legend.position = "bottom")
ggsave("./figures/fig5_pearson.png")


