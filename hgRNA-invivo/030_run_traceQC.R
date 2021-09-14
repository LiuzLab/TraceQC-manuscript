library(readr)
library(dplyr)
library(stringr)
library(tools)
library(TraceQC)
library(fastqcr)
library(readr)

sra <- read_csv("./data/000_SraRunTable.txt") %>%
  select(Run,`Library Name`)

qc_dir <- "./fastqc/"
fastq_dir <- "./data/020_fastq_by_identifier"
for (dir in list.dirs(fastq_dir)){
    fastqc(dir,qc.dir=qc_dir) 
}

identifiers <- read_csv("./data/000_ref/hgRNA_identifiers.csv")
all_length <- c(21,25,30,35)
for (l in all_length) {
  ref <- read_lines(sprintf("./data/000_ref/L%s.txt",l))
  refseq <- ref[1]
  regions_str <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions_str) %>%
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(.data$start),
           end=strtoi(.data$end)) %>%
    mutate(region=as.character(.data$region))
  # tmp <- list(refseq=refseq,regions=regions)
  # plot_construct(tmp)
  
  L_identifiers <- filter(identifiers,Length>=(l-1)&Length<=(l+1))
  for (i in 1:nrow(L_identifiers)) {
    identifier <- as.character(L_identifiers[i,"Identifier (ID)"])
    spacer <- as.character(L_identifiers[i,"Spacer regions (TSS to PAM)"])
    spacer_start <- regions[regions$region=="spacer","start"]
    spacer_end <- regions[regions$region=="spacer","end"]
    ref_id_seq <- paste(substr(refseq,start=1,stop=spacer_start-1),spacer,
                        substr(refseq,start=spacer_end+1,stop=nchar(refseq)),sep="")
    # tmp <- list(refseq=ref_id_seq,regions=regions)
    # plot_construct(tmp)
    out_file <- sprintf("./data/000_ref/L%s_%s.txt",l,identifier)
    write(paste(c(ref_id_seq,ref[2:3]),sep="\n"),
          out_file)}}


for (f in list.files(fastq_dir,recursive=TRUE)) {
  tmp <- strsplit(f,"/")[[1]][2]
  tmp <- strsplit(file_path_sans_ext(tmp),split="_")[[1]]
  sra <- tmp[1]
  identifier <- tmp[2]
  
  input_file <- sprintf("%s/%s",fastq_dir,f)
  ref_file <- list.files("./data/000_ref/",pattern=identifier)
  ref_file <- sprintf("./data/000_ref/%s",ref_file)
  output_file <- sprintf("./data/030.1_alignment/%s_%s.txt",sra,identifier)
  sequence_alignment(input_file=input_file,ref_file=ref_file,
                     output_file=output_file)}

for (input_file in list.files("./data/030.1_alignment/")) {
  print(input_file)
  aligned_reads <- read_tsv(sprintf("./data/030.1_alignment/%s",input_file))
  if (nrow(aligned_reads)>0) {
  traceQC_input <- list(aligned_reads=aligned_reads)
  mutation_event <- seq_to_character(traceQC_input,ncores=1,
                                     use_CPM=FALSE,alignment_score_cutoff=-Inf,
                                     abundance_cutoff=0)
  write_tsv(mutation_event,sprintf("./data/030.2_traceQC_obj/%s",input_file))}
}

for (ref_file in list.files("./data/000_ref/",pattern="L[0-9][0-9]_")) {
  alignment_threshold <- sequence_permutation(ref_file=sprintf("./data/000_ref/%s",ref_file))
  write_tsv(alignment_threshold,sprintf("./data/030.3_alignment_threshold/alignment_threshold_%s",ref_file))}
