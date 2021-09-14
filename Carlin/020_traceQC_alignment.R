library(TraceQC)
library(fastqcr)
library(tools)
library(stringr)
library(readr)
library(ggpubr)
library(dplyr)

ref_file <- "./data/000_ref.txt"
fastq_dir <- "./data/011_merged_fastq/"

for (f in list.files(fastq_dir)) {
  name <- tools::file_path_sans_ext(f)
  out_alignment_file <- paste("./data/020_alignment/",name,".txt",sep="")
  if (!file.exists(out_alignment_file)) {
    print(name)
    sequence_alignment(input_file=paste(fastq_dir,f,sep=""),
                       ref_file=ref_file,output_file=out_alignment_file,ncores=16)}}

df <- sequence_permutation(ref_file=ref_file)
write_tsv(df,"./data/021_alignment_threshold.txt")

library(reticulate)
get_abspath <- function(f) {
  if(is.null(f)) {
    stop("NULL shouldn't be an argument.")
  }
  file.path(normalizePath(dirname(f)), basename(f))
}
input_file="./data/tmp.fastq"
ref_file <- "./data/000_ref.txt"
output_file="./tmp.txt"
ncores=1
match=2
mismatch=-2
gapopen=-6
gapextension=-0.1
ncores = 4
args <- list("input"=get_abspath(input_file),
               "reference"=get_abspath(ref_file),
               "output"=get_abspath(output_file),
               "match"=match,"mismatch"=mismatch,"gapopen"=gapopen,
               "gapextension"=gapextension,
               "ncores"=ncores)
source_python(system.file("py", "alignment.py", package="TraceQC"))
module <- reticulate::import_from_path("alignment", system.file("py", package = "TraceQC"))
  
  
tic("Alignment")
message(paste0("Running an alignment between ", ref_file,
               " and ", input_file, "."))
module$alignment(args)
toc()
  
