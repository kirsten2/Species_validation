#!/usr/bin/env python3
import sys
import glob
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

#Usage: python3 concatenate_aln.py 

#Identify the files that are to be concatenated (using file suffix *aln_prune.fasta. Check that files have been identified.

aln_files = glob.glob('./*aln_trim.fasta')
if (len(aln_files) == 0):
    print('No alignment files with suffix "*aln_prune.fasta" detected, nothing to concatenate')
    print('Exiting script')
    exit()
print('Concatenating ',len(aln_files), 'alignment files')

#Initiate new alignment object to hold the concatenate alignment. Sort each alignment before appending
new_aln = None
for file in aln_files:
    aln = AlignIO.read(file, "fasta")
    aln.sort()
    if (new_aln == None):
        new_aln = aln
    else:
        new_aln = new_aln + aln

#Generate the output-file
fh_out = open('cat_all.phy',"w")
AlignIO.write(new_aln,fh_out,"phylip")

