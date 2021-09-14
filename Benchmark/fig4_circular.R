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
    max_cnt <- ceiling(max(df$log10_count))
    for (i in 1:nrow(ref$regions)) {
      if (ref$regions$region[i] %in% names(annotation_cols)) {
        circos.rect(xleft=ref$regions$start[i]-target_start,ybottom=0,
                    xright=ref$regions$end[i]-target_start+1,ytop=1,
                    col=annotation_cols[ref$regions$region[i]],border=NA)}}
    
    col_intra <- colorRamp2(c(0,ceiling(max(df$log10_count))), 
                                c("white", "blue2"))
    col_inter <- colorRamp2(c(0,ceiling(max(df$log10_count))), 
                                c("white", "red2"))
    
    df <- mutate(df,color=case_when(deletion_type=="intra-target" ~ col_intra(log10_count),
                                    deletion_type=="inter-target" ~ col_inter(log10_count)))
    for (i in 1:nrow(df)) {
    circos.link(1,df$start[i],1,df$end[i],
                  h.ratio = 0.9,lwd = 1,
                  col = alpha(df$color[i], df$log10_count[i] / max_cnt))}}

################ Carlin ################ 
Carlin <- readRDS("../Carlin/data/030_all_data.rds")
ref <- parse_ref_file("../Carlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
ref_link_end <- c(ref_link_end)-ref_target_start+2

deletion <- Carlin$ChronicInduction_NegDox_Trial1_0h %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end)+1,
         ending_sgRNA=findInterval(start+length,ref_link_end)+1) %>%
  mutate(deletion_type=case_when(starting_sgRNA==ending_sgRNA ~ "intra-target",
                                 TRUE ~ "inter-target")) %>%
  mutate(starting_cutsite=start+1-ref_link_end[starting_sgRNA+1],
         ending_cutsite=start+length+1-ref_link_end[ending_sgRNA+1],
         end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_carlin.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()

################ scCarlin ################
scCarlin <- readRDS("../scCarlin/data/050.1_all_data.rds")
ref <- parse_ref_file("../scCarlin/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
ref_link_end <- c(ref_link_end)-ref_target_start+2

deletion <- scCarlin$`5FU_FO837_SC_Amplicon` %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=n()) %>%
  ungroup %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end)+1,
         ending_sgRNA=findInterval(start+length,ref_link_end)+1) %>%
  mutate(deletion_type=case_when(starting_sgRNA==ending_sgRNA ~ "intra-target",
                                 TRUE ~ "inter-target")) %>%
  mutate(starting_cutsite=start+1-ref_link_end[starting_sgRNA+1],
         ending_cutsite=start+length+1-ref_link_end[ending_sgRNA+1],
         end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_scCarlin.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()

################ GESTALT ################
GESTALT <- readRDS("../GESTALT/data/030.1_all_data_whole_organism.rds")
ref <- parse_ref_file("../GESTALT/data/000_ref_v6.txt")
ref_link_end <- filter(ref$regions,region=="Link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
ref_link_end <- c(ref_link_end)-ref_target_start+2

deletion <- GESTALT$SRR3561150 %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end)+1,
         ending_sgRNA=findInterval(start+length,ref_link_end)+1) %>%
  mutate(deletion_type=case_when(starting_sgRNA==ending_sgRNA ~ "intra-target",
                                 TRUE ~ "inter-target")) %>%
  mutate(starting_cutsite=start+1-ref_link_end[starting_sgRNA+1],
         ending_cutsite=start+length+1-ref_link_end[ending_sgRNA+1],
         end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_GESTALT.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()

################ scGESTALT ################
scGESTALT <- readRDS("../scGESTALT/data/050.1_all_data.rds")
ref <- parse_ref_file("../scGESTALT/data/000_ref.txt")
ref_link_end <- filter(ref$regions,region=="link") %>%
  pull(end)
ref_target_start <- filter(ref$regions,region=="target") %>%
  pull(start)
ref_target_end <- filter(ref$regions,region=="target") %>%
  pull(end)
ref_link_end <- c(ref_link_end)-ref_target_start+2

deletion <- scGESTALT$SRR6176748 %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=n()) %>%
  ungroup %>%
  mutate(starting_sgRNA=findInterval(start,ref_link_end)+1,
         ending_sgRNA=findInterval(start+length,ref_link_end)+1) %>%
  mutate(deletion_type=case_when(starting_sgRNA==ending_sgRNA ~ "intra-target",
                                 TRUE ~ "inter-target")) %>%
  mutate(starting_cutsite=start+1-ref_link_end[starting_sgRNA+1],
         ending_cutsite=start+length+1-ref_link_end[ending_sgRNA+1],
         end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_scGESTALT.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()

################ hgRNA-invitro ################
hgRNA_invitro <- readRDS("../hgRNA-invitro/data/040.1_all_data.rds")
ref <- parse_ref_file("../hgRNA-invitro/data/000_ref/A21.txt")
deletion <- hgRNA_invitro$`A21-0d(A')` %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  mutate(deletion_type="intra-target") %>%
  mutate(end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_hgRNA-invitro.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()

################ hgRNA-invivo ################
hgRNA_invivo <- readRDS("../hgRNA-invivo/data/040.1_all_data.rds")
ref <- parse_ref_file("../hgRNA-invivo/data/000_ref/L21.txt")
identifiers <- c()
for (ref_file in list.files("../hgRNA-invivo/data/000_ref/",pattern="L21_")) {
  identifier <- file_path_sans_ext(ref_file)
  identifiers <- c(identifiers,strsplit(identifier,split="_")[[1]][2])}

deletion <- hgRNA_invivo$SRR7633623 %>%
  filter(type=="deletion",identifier %in% identifiers) %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=sum(count)) %>%
  ungroup %>%
  mutate(deletion_type="intra-target") %>%
  mutate(end=start+length,log10_count=log10(count+1))

png("./figures/fig4_circular_hgRNA_invivo.png",width=300,height=300)
circular_chordgram(deletion,ref)
dev.off()
