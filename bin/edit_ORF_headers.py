#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Usage: edit_ORF_headers.py sample.fasta

#Parse sample-name from input file-name, and generate outfile-name
infile_name = sys.argv[1]
infile_name_split = infile_name.split('.')
sample = infile_name_split[0]
outfile = sys.stdout #Writing to stdout 

#Try open the input-file, or exit
try:
    infile = sys.argv[1]
    fh_infile = open(infile)
except:
    print('The following fasta-file couldnt be opened:')
    print(infile)
    print('Exiting script')
    exit()

#Generate new seq-records with sample-name added to the header
new_seq_records = list()
for record in SeqIO.parse(fh_infile, "fasta"):
    old_header = record.id
    new_header = sample + '_' + old_header
    seq = record.seq
    new_record = SeqRecord(seq, id=new_header, description='')
    new_seq_records.append(new_record)

#Print new seq-records to fasta-file
SeqIO.write(new_seq_records, outfile, "fasta")