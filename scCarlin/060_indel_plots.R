library(readr)
library(dplyr)
library(tidyr)
library(TraceQC)
library(ggplot2)

deletion <- read_tsv("./data/040_cell_mutations/5FU_FO817_SC_Amplicon.txt") %>%
  filter(type=="deletion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=n()) %>%
  ungroup
ref <- parse_ref_file("./data/000_ref.txt")
plot_construct(ref)
ggsave("./figures/supfig2_v1_construct.png")

traceQC_input <- list(refseq=ref$refseq,regions=ref$regions,mutation=deletion)
png("./figures/fig3_deletion_pattern.png")
plot_deletion_hotspot(traceQC_input)
dev.off()

insertion <- read_tsv("./data/040_cell_mutations/5FU_FO817_SC_Amplicon.txt") %>%
  filter(type=="insertion") %>%
  group_by(type,start,length,mutate_to) %>%
  summarise(count=n()) %>%
  ungroup
traceQC_input <- list(refseq=ref$refseq,regions=ref$regions,mutation=insertion)
png("./figures/fig3_insertion_pattern.png")
plot_insertion_hotspot(traceQC_input)
dev.off()
