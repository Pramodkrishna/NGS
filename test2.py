#!/usr/bin/python
import Bio
from Bio import SeqIO

"""
for seq_record in SeqIO.parse("sample_two.fastq", "fastq"):
    #rec_seq  = seq_record.id
    rec_seq  = repr(seq_record.seq)
    phred_score =  seq_record.letter_annotations["phred_quality"]
    print(phred_score)
    print(rec_seq)
"""
#sub_rec = seq_record[5:15]
for seq_record in SeqIO.parse("fastaQ.fastq", "fastq"):
    #rec_seq  = seq_record.id
    #rec_seq  = repr(seq_record.seq)
    #phred_score =  seq_record.letter_annotations["phred_quality"]
    sub_rec = seq_record[5:15]
    phred_score = sub_rec.letter_annotations["phred_quality"]
    print(phred_score)
   #print(seq_record.format("qual"))
    	
   # print(rec_seq)
	
