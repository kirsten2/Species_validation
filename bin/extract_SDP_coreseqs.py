#!/usr/bin/env python3
import sys
import os
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

#Usage: python3 bin/extract_SDP_coreseqs.py Candidate_SDPs.txt firm5 

def get_SDP_genome(fh_file):
    SDP_genome_dict = dict()
    for line in fh_file:
        if (line.find('SDP') != -1): continue
        line = line.strip()
        split_line = line.split('\t')
        genome_id = split_line[0]
        SDP = split_line[2]
        if (SDP not in SDP_genome_dict):
            SDP_genome_dict[SDP]=dict()
        SDP_genome_dict[SDP][genome_id]=1
    fh_file.close()
    return(SDP_genome_dict)

def get_genome_id_seq_obj(filename):
    seq_record_dict = dict()
    for seq_record in SeqIO.parse(filename, "fasta"):
        gene_id = seq_record.id
        split_id = gene_id.split('_')
        genome_id = split_id[0]
        seq_record_dict[genome_id]=seq_record
    return(seq_record_dict)

#Read input-file with genome-ids of candidate SDPs for validation. Get SDPs from the file, and prepare names of SDP output directories
try:
    cand_SDP = sys.argv[1]
    fh_cand_SDP = open(cand_SDP)
except:
    print('Provide tab-delimited text-file with candidate SDPs. Exiting script')
    exit()

SDP_genome = get_SDP_genome(fh_cand_SDP)
SDPs = list(SDP_genome.keys())
SDPs.sort()
current_dir = os.getcwd()
SDP_dirs = dict()
for SDP in SDPs:
    SDP_dir = current_dir + '/' + SDP
    SDP_dirs[SDP] = SDP_dir

#Check that the user-provided directory for alignments that will be subset exists
dir_in = sys.argv[2]
isdir = os.path.isdir(dir_in)
if (isdir == False):
    print('The directory provided for the input files doest exist, check the name/path:')
    print(dir_in)
    print('Exiting script!')
    exit()

#Move to the input alignment directory, get the names of files to be subset (or exit, if there are no such files)
os.chdir(dir_in)
ffn_suffix = '*ffn'
ffn_files = glob.glob(ffn_suffix)
nb_ffn_files = len(ffn_files)
if (nb_ffn_files == 0):
     print('No fasta-files with the suffix "*ffn" were found in the SDP directory. Check the directory:')
     print(dir_in)
     print('Exiting script!')
     exit()

#Create the SDP output directories (or exit, if they already exist)
for SDP in SDP_dirs.keys():
    try: 
        os.mkdir(SDP_dirs[SDP])
    except OSError:
        print('Could not create the SDP output directory:') 
        print(SDP_dirs[SDP])
        print('Exiting script!')   
        exit()

#Subset the alignment files by SDP, and print to the corresponding output directories
count_file_progress = 0
for ffn_file in ffn_files:
    seq_records = get_genome_id_seq_obj(ffn_file)
    count_file_progress += 1    
    for SDP in SDP_genome.keys():
        genomes = list(SDP_genome[SDP].keys())
        SDP_seq_records = {key: seq_records[key] for key in genomes}
        SDP_ffn_outfile = SDP_dirs[SDP] + '/' + ffn_file
        SeqIO.write(SDP_seq_records.values(), SDP_ffn_outfile, "fasta")
    if (count_file_progress % 100 == 0):
        print('Done sub-setting', count_file_progress, ' ffn-files..')

#Check the content of the SDP directories, print a file with SDP-dirs containing subset fasta files. Print information on the SDP-directories created to stdout. 

os.chdir(current_dir)
final_SDP_dirs = list()
for SDP in SDP_dirs.keys():
    os.chdir(SDP_dirs[SDP])
    ffn_files_SDP = glob.glob(ffn_suffix)
    nb_ffn_files_SDP= len(ffn_files)
    if (nb_ffn_files != 0):
        final_SDP_dirs.append(SDP)
    os.chdir(current_dir)

if (len(final_SDP_dirs) > 0):
    final_SDP_dirs.sort()
    fh_out = open('SDP_dirs.txt','w')
    print('The following SDP-directories with subset data were succesfully created:')
    for SDP in final_SDP_dirs:
        print(SDP)
        fh_out.write(SDP_dirs[SDP] + '\n')
    fh_out.close()
else:
    print('ERROR: All the created SDP directories are empty! Check that the alignment files contain SDPs with genome-ids corresponding precisely the data in "Candidate_SDPs.txt" file')
    
        


