#!/usr/bin/env python3
import sys
import os
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

#Usage: python3 bin/extract_SDP_coreseqs.py Candidate_species.txt firm5 

def get_species_genome(fh_file):
    species_genome_dict = dict()
    line_count = 0
    for line in fh_file:
        line_count = line_count + 1
        if (line_count == 1): continue      
        line = line.strip()
        split_line = line.split('\t')
        genome_id = split_line[0]
        species = split_line[2]
        if (species not in species_genome_dict):
            species_genome_dict[species]=dict()
        species_genome_dict[species][genome_id]=1
    fh_file.close()
    return(species_genome_dict)

def get_genome_id_seq_obj(filename):
    seq_record_dict = dict()
    for seq_record in SeqIO.parse(filename, "fasta"):
        gene_id = seq_record.id
        split_id = gene_id.split('_')
        genome_id = split_id[0]
        seq_record_dict[genome_id]=seq_record
    return(seq_record_dict)

#Read input-file with genome-ids of candidate SDPs for validation. Get SDPs from the file, and prepare names of SDP output directories
fh_log_out = open('log.txt','a')
fh_cand_species = open(sys.argv[1])
species_genome = get_species_genome(fh_cand_species)
species = list(species_genome.keys())
species.sort()
run_dir = os.getcwd()
species_dirs = dict()
for spec in species:
    species_dir = run_dir + '/' + spec
    species_dirs[spec] = species_dir

#Check that the user-provided directory for alignments that will be subset exists
dir_in = sys.argv[2]
isdir = os.path.isdir(dir_in)
if (isdir == False):
    fh_log_out.write('The directory provided for the input files doest exist, check the name/path:')
    fh_log_out.write(dir_in + "\n")
    exit()

#Move to the input alignment directory, get the names of files to be subset (or exit, if there are no such files)
os.chdir(dir_in)
ffn_suffix = '*ffn'
ffn_files = glob.glob(ffn_suffix)
nb_ffn_files = len(ffn_files)
if (nb_ffn_files == 0):
    fh_log_out.write('No fasta-files with the suffix "*ffn" were found in the input directory. Check the directory:\n')
    fh_log_out.write(dir_in + "\n")
    exit()

#Create the SDP output directories
for spec in species_dirs.keys():
    try: 
        os.mkdir(species_dirs[spec])
    except OSError:
        fh_log_out.write('Could not create the species output directory:\n') 
        fh_log_out.write(species_dirs[spec] + "\n")

#Subset the alignment files by SDP, and print to the corresponding output directories
count_file_progress = 0
for ffn_file in ffn_files:
    seq_records = get_genome_id_seq_obj(ffn_file)
    count_file_progress += 1    
    for spec in species_genome.keys():
        genomes = list(species_genome[spec].keys())
        species_seq_records = {key: seq_records[key] for key in genomes}
        species_ffn_outfile = species_dirs[spec] + '/' + ffn_file
        if (os.path.isfile(species_ffn_outfile)):
            fh_log_out.write('The following subset ffn-file already exists:' +  species_ffn_outfile + "\n")
        else:
            SeqIO.write(species_seq_records.values(), species_ffn_outfile, "fasta")
    if (count_file_progress % 100 == 0):
        print('Done sub-setting', count_file_progress, ' ffn-files..')

#Return to run directory, and print a few stats to summary file
os.chdir(run_dir)
fh_log_out.write("Successfully created the following directories:" + "\n")
for spec in species_dirs.keys():
    fh_log_out.write(spec + "\n")
fh_log_out.write("A total of " + str(nb_ffn_files) + " core gene family ffn-files were subset to each directory\n\n")
fh_log_out.close()

