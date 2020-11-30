#!/usr/bin/env python3
import sys
import os
import glob
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq

#Usage: python3 recruit_ORFs.py ORFs_db_example.ffn SDP_dir

def get_hit_ids(blast_file):
    blastfile_handle = open(blast_file)
    blast_records = NCBIXML.parse(blastfile_handle)
    hit_ids = list()
    for blast_record in blast_records:
        for alignment in blast_record.alignments:              
            hit_name_full = alignment.title              
            hit_name_split = hit_name_full.split()              
            hit_id=hit_name_split[-1]
            hit_ids.append(hit_id)
    return(hit_ids)

#Parse input arguments, set of working dirs
orf_db_file = sys.argv[1]
sdp_dir = sys.argv[2]
current_dir = os.getcwd()
fh_log_out = open('log.txt','a')

#Move into SDP directory, use glob to get the names of all blast-files (exit script, if there are no such files)
os.chdir(sdp_dir)
blast_suffix = '*blastn'
blast_files = glob.glob(blast_suffix)
nb_blast_files = len(blast_files)
if (nb_blast_files == 0):
    print('ERROR: No blast-files with the suffix "*blastn" were found in the SDP directory: \n')
    fh_log_out.write('ERROR: No blast-files with the suffix "*blastn" were found in the SDP directory: \n')
    fh_log_out.write(sdp_dir + "\n")
    exit()

#Loop over blast-files, and get the hit-ids for all hits. Store in dictionaries, with OG-affiliation. 
hit_id_OG = dict()
OG_hit_id = dict()
count_progress = 0
for file in blast_files:
    count_progress += 1
    if (os.stat(file).st_size != 0): #check if file is empty
        split_filename = file.split('.')
        OG=split_filename[0]
        hit_ids = get_hit_ids(file)
        OG_hit_id[OG]=dict()
        for id in hit_ids:
            hit_id_OG[id]=OG
            OG_hit_id[OG][id]=1
    if (count_progress % 100 == 0):
        print('Finished parsing',count_progress,'blast-files')

#Move back to run-directory, parse out seq-objects for all hit-ids
os.chdir(current_dir)
hit_seq_objects = dict()
print('Getting seq-records from ORF db')
for seq_record in SeqIO.parse(orf_db_file, "fasta"):
    if(seq_record.id in hit_id_OG):
        seq_length = len(seq_record.seq)
        if (seq_length > 200):
            hit_seq_objects[seq_record.id] = seq_record

            
#Print ORF sequences to the SDP directory
os.chdir(sdp_dir)
print('Printing recruited ORFs to files')
nb_orfs_recruited=0
for OG in OG_hit_id.keys():
    ORF_seq_objects = list()
    for hit_id in OG_hit_id[OG].keys():
        if (hit_id in hit_seq_objects):
            ORF_seq_objects.append(hit_seq_objects[hit_id])
    nb_objects = len(ORF_seq_objects)
    if (nb_objects != 0):
       orf_outfile = OG + '_orfs.ffn'
       if (os.path.isfile(orf_outfile)):
           fh_log_out = open('log.txt','a')
           fh_log_out.write("NOTE: The following orf-file already exists: " + orf_outfile + "\n")
           fh_log_out.close()
           print('NOTE: The following orf-file already exists:', orf_outfile)
       else:
           nb_orfs_recruited += len(ORF_seq_objects)
           SeqIO.write(ORF_seq_objects, orf_outfile, "fasta")
os.chdir(current_dir)
fh_log_out = open('log.txt','a')
fh_log_out.write("Recruited " + str(nb_orfs_recruited) + " ORFs to dir: " + sdp_dir + "\n")
fh_log_out.close()
print("Recruited",nb_orfs_recruited, "ORFs to SDP", sdp_dir)
