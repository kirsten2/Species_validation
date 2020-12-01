#!/usr/bin/env python3
import sys
import os
import glob

#Usage: python3 get_species_dirs.py Candidate_species.txt

#Read Candidate species file
try:
    open(sys.argv[1])
except:
    print('Cant open candidate species file:')
    print(sys.argv[1])
    print('Exiting script!')
    exit()

cand_species_dict = dict()
fh_cand_species = open(sys.argv[1])
for line in fh_cand_species:
    line=line.strip()
    split_line = line.split('\t')
    cand_species_dict[split_line[2]] = 1
fh_cand_species.close()

#Check that each candidate species directory was successfully created and contains fasta-files for processing

species_dirs = cand_species_dict.keys()
species_dirs_verified = list()
run_dir = os.getcwd()
for dir in species_dirs:
    isdir = os.path.isdir(dir)
    if (isdir == True):
        os.chdir(dir)
        ffn_suffix = '*ffn'
        ffn_files = glob.glob(ffn_suffix)
        try:
            open(ffn_files[0])
            species_dirs_verified.append(dir)
        except Exception:
            pass
        os.chdir(run_dir)

#Print verified species-directories as string to stdout for further processing with bash script

if (len(species_dirs_verified) > 0):
    species_dirs_verified.sort()
    species_dir_string = " ".join(species_dirs_verified)
    print(species_dir_string)
else:
    print('ERROR')