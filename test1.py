#!/usr/bin/python
import Bio
from Bio import SeqIO

"""
fasta_sequences = SeqIO.parse(open("sample_one.fastaq"),'fastaq')
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)
#parse a fastaq file and open it

"""

"""
fasta_file = open("sample_one.fastaq")
fas_read = fasta_file.read()
print fas_read
#Read a fasta file without using Bio module
"""

"""
for record in SeqIO.parse("sample_one.fastq", "fastq"):
	 #print("%s %s" % (record.id, record.seq))
	 #print(record)
	 print(record.letter_annotations["phred_quality"])

"""

sqe_file = open("fasta.txt")
cont = sqe_file.read()
#print cont
#for line in cont:
    #positions = cont.split()
    #positions_one = positions[0]
    #positions_two = positions[1]
    #combined = positions_one + positions_two
    #print positions
    #print combined
cont_line = cont.split()
#print cont_line
contstr1 = cont_line[0]
contstr2 = cont_line[1]
#combined =  contstr1 + contstr2
#print combined
constr_complete = cont_line[10:]
print	' , '.join(constr_complete)

#String or data manipulation 
#Open the fasta file, split it into different lines which returns a list 
#Access each element in the list
#Join them at the end and return a complete list of elements



