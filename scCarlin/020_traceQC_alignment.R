library(TraceQC)

ref_file <- "./data/000_ref.txt"
ref <- parse_ref_file(ref_file)
plot_construct(ref)
for (name in list.files("./data/010_cellranger_output/",pattern="Amplicon")) {
  input_file <- sprintf("./data/010_cellranger_output/%s/outs/possorted_genome_bam.bam",name)
  out_alignment_file <- paste("./data/020_alignment_10x/",name,".txt",sep="")
  sequence_alignment_for_10x(input_file,ref_file,output_file=out_alignment_file,ncores=10)}
