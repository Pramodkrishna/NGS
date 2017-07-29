#!/usr/bin/python
import Bio
from Bio import SeqIO
#def phred_score(A):
for record in SeqIO.parse("SRR494011_1.fastq", "fastq"):
    #print("%s %s" % (record.id, record.seq))
    #print(record)
    #print(record.letter_annotations["phred_quality"])

