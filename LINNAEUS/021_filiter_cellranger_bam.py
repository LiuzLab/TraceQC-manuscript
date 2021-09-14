import pysam
from glob import glob

if __name__ == "__main__":
    for input_file in glob("./data/020_cellranger_output/*/outs/possorted_genome_bam.bam"):
        run = input_file.split("/")[3]
        print(run)
        output_file = "./data/021_cellranger_bam_filtered/%s.bam"% run
        inputbam = pysam.AlignmentFile(input_file, "rb")
        filteredreads = pysam.AlignmentFile(output_file, "wb", template=inputbam)
        for read in inputbam.fetch(until_eof=True):
            if read.is_unmapped: 
                try:
                    read.get_tag("CB")
                except:
                    continue
                try:
                    read.get_tag("UB")
                except:
                    continue
                filteredreads.write(read)

        filteredreads.close()
        inputbam.close()
