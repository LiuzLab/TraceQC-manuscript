import os
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from multiprocessing import Process, Pool
from tqdm import tqdm
import gzip
import subprocess

class align_identifier(object):
    def __init__(self,srr):
        identifiers = pd.read_csv("./000_ref/hgRNA_identifiers.csv")
        identifiers = list(identifiers["Identifier (ID)"])
        self.srr = srr
        self.r1 = SeqIO.parse("./010_raw/%s_1.fastq"%srr,"fastq")
        self.r2 = SeqIO.parse("./010_raw/%s_2.fastq"%srr,"fastq")
        self.signature_sequence = "GGGCCCGAATTC"
        self.identifiers = identifiers
        self.identifiers = dict(zip(identifiers,[[] for i in identifiers]))
    
    def get_identifier(self,param):
        r1,r2 = param
        aligned_seq = pairwise2.align.localms(r2.seq,self.signature_sequence,1,-1,-1,-1,one_alignment_only=True)[0]
        if aligned_seq.score >= 10:
            identifier = r2.seq._data[aligned_seq.end:aligned_seq.end+10]
            if identifier in self.identifiers:
                self.identifiers[identifier].append(r1)

    def alignment(self):
        res = []
        for r in zip(self.r1,self.r2):
            self.get_identifier(r)
    
    def write_files(self):
        subprocess.call(["mkdir","-p","./020_fastq_by_identifier/%s"%self.srr])
        
        for identifier in self.identifiers:
            if len(self.identifiers[identifier]) > 0:
                with open("./020_fastq_by_identifier/%s/%s_%s.fastq"%(self.srr,self.srr,identifier), "w") as output_handle:
                    SeqIO.write(self.identifiers[identifier], output_handle, "fastq")
                    
def align_id(srr):
    myClass = align_identifier(srr)
    myClass.alignment()
    myClass.write_files()
    
sra = pd.read_csv("./000_SraRunTable.txt")
runs = list(sra["Run"])

num_threads = 16
pool = Pool(num_threads)
pool.map(align_id,runs)
