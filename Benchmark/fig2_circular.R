library(TraceQC)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)

circular_chordgram <-
  function(df, ref) {
    regions <- ref$regions
    target_start <- regions %>% filter(.data$region == "target") %>% pull(.data$start)
    target_end <- regions %>% filter(.data$region == "target") %>% pull(.data$end)
    refseq <- substr(ref$refseq, start = target_start, stop = target_end)
    l <- nchar(refseq)+1
    
    circos.par(start.degree = 90)
    circos.initialize(factors = factor(1), xlim = c(0, ceiling(l * 1.02)))
    circos.track(track.height = 0.15, ylim = c(0, 1), bg.border = NA)
    circos.rect(xleft=0,ybottom=0,xright=l+1,ytop=1,col="#E5E7E9",border=NA)
    annotation_cols <- c("sgRNA"="#CCD1D1","pam"="#17202A","spacer"="#CCD1D1","Pam"="#17202A")
    for (i in 1:nrow(ref$regions)) {
      if (ref$regions$region[i] %in% names(annotation_cols)) {
        circos.rect(xleft=ref$regions$start[i]-target_start,ybottom=0,
                    xright=ref$regions$end[i]-target_start+1,ytop=1,
                    col=annotation_cols[ref$regions$region[i]],border=NA)}}
    
    indel_colors <- c("deletion"="#F8766D","insertion"="#00BFC4","substitution"="#7CAE00")
    for (i in 1:nrow(df)) {
      if (df$type[i]=="substitution") {
        circos.rect(xleft=df$start[i],ybottom=0,
                    xright=df$end[i],ytop=1,
                    col=indel_colors[df$type[i]],border=NA)
        } else {circos.link(1,df$start[i],1,df$end[i],
            h.ratio = 0.9,lwd = 2,col = indel_colors[df$type[i]])}}}

################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
indels <- data.frame()
for (sample_name in names(Carlin)[2:length(Carlin)]) {
  tmp <- filter(Carlin[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=sum(count)) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
    indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
      mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
    indels$count <- indels$count.x + indels$count.y
    indels$count.x <- NULL; indels$count.y <- NULL
    }
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

png("./figures/fig2_circular_carlin.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
ref <- parse_ref_file("../scCarlin/data/000_ref.txt")
indels <- data.frame()
for (sample_name in names(scCarlin)[2:length(scCarlin)]) {
  tmp <- filter(scCarlin[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=n()) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
      indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
        mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
      indels$count <- indels$count.x + indels$count.y
      indels$count.x <- NULL; indels$count.y <- NULL
    }
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length)
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length)

png("./figures/fig2_circular_sccarlin.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
indels <- data.frame()
for (sample_name in names(GESTALT)[2:length(GESTALT)]) {
  tmp <- filter(GESTALT[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=sum(count)) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
      indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
        mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
      indels$count <- indels$count.x + indels$count.y
      indels$count.x <- NULL; indels$count.y <- NULL
    }
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

png("./figures/fig2_circular_GESTALT.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
ref <- parse_ref_file("../scGESTALT/data/000_ref.txt")
indels <- data.frame()
for (sample_name in names(scGESTALT)[2:length(scGESTALT)]) {
  tmp <- filter(scGESTALT[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=n()) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
      indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
        mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
      indels$count <- indels$count.x + indels$count.y
      indels$count.x <- NULL; indels$count.y <- NULL
    }
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length)
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length)

png("./figures/fig2_circular_scGESTALT.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ hgRNA-invitro ################
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
ref <- parse_ref_file("../hgRNA-invitro/data/000_ref/A21.txt")
indels <- data.frame()
for (sample_name in names(hgRNA_invitro)[2:length(hgRNA_invitro)]) {
  if (str_detect(sample_name,"A21")) {
  tmp <- filter(hgRNA_invitro[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=sum(count)) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
      indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
        mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
      indels$count <- indels$count.x + indels$count.y
      indels$count.x <- NULL; indels$count.y <- NULL
    }}
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

png("./figures/fig2_circular_hgRNA-invitro.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
ref <- parse_ref_file("../hgRNA-invivo/data/000_ref/L21.txt")
identifiers <- c()
for (ref_file in list.files("../hgRNA-invivo/data/000_ref/",pattern="L21_")) {
  identifier <- file_path_sans_ext(ref_file)
  identifiers <- c(identifiers,strsplit(identifier,split="_")[[1]][2])}

indels <- data.frame()
for (sample_name in names(hgRNA_invivo)[2:length(hgRNA_invivo)]) {
  if (nrow(hgRNA_invivo[[sample_name]])>100) {
    tmp <- filter(hgRNA_invivo[[sample_name]],type != "unmutated") %>%
      filter(identifier %in% identifiers) %>%
      group_by(type,start,length) %>%
      summarise(count=sum(count)) %>%
      ungroup
    if (nrow(indels)==0) {
      indels <- tmp} else {
        indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
          mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
        indels$count <- indels$count.x + indels$count.y
        indels$count.x <- NULL; indels$count.y <- NULL
      }}
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

png("./figures/fig2_circular_hgRNA_invivo.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()

################ LINNAEUS ################
LINNAEUS <- readRDS("../LINNAEUS/data/060.1_all_data.rds")
ref <- parse_ref_file("../LINNAEUS/data/000_ref.txt")
target_start <- filter(ref$regions,region=="target") %>% pull(start)
target_end <- filter(ref$regions,region=="target") %>% pull(end)
sgRNA_start <- filter(ref$regions,region=="sgRNA") %>% pull(start)
sgRNA_end <- filter(ref$regions,region=="sgRNA") %>% pull(end)
ref$regions <- mutate(ref$regions,region=case_when(region=="sgRNA" ~ "spacer",
                                                   TRUE ~ region))
ref$regions$start <- c(sgRNA_start-20,sgRNA_start,sgRNA_start+17)


indels <- data.frame()
for (sample_name in names(LINNAEUS)[2:length(LINNAEUS)]) {
  tmp <- filter(LINNAEUS[[sample_name]],type != "unmutated") %>%
    group_by(type,start,length) %>%
    summarise(count=n()) %>%
    ungroup
  if (nrow(indels)==0) {
    indels <- tmp} else {
      indels <- full_join(indels,tmp,by=c("type","start","length")) %>%
        mutate(count.x=replace_na(count.x,0),count.y=replace_na(count.y,0))
      indels$count <- indels$count.x + indels$count.y
      indels$count.x <- NULL; indels$count.y <- NULL
    }
}
df <- arrange(indels,desc(count)) %>%
  head(30) %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))
df <- group_by(indels,type) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  ungroup %>%
  mutate(end=start+length) %>%
  mutate(type=case_when(type=="mutation" ~ "substitution",
                        TRUE ~ type))

png("./figures/fig2_circular_LINNAEUS.png",width=300,height=300)
circular_chordgram(df,ref)
dev.off()


