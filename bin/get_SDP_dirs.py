#!/usr/bin/env python3
import sys
import os
import glob

#Usage: python3 get_SDP_dirs.py Candidate_SDPs.txt

#Read Candidate SDP file
try:
    open(sys.argv[1])
except:
    print('Cant open candidate SDP-file:')
    print(sys.argv[1])
    print('Exiting script!')
    exit()

cand_SDP_dict = dict()
fh_cand_SDP = open(sys.argv[1])
for line in fh_cand_SDP:
    line=line.strip()
    split_line = line.split('\t')
    try:
        cand_SDP_dict[split_line[2]] = 1
    except:
        print('ERROR Candidate SDP-file appear to be wrongly formatted, there is no data in the third column, exiting bash-script')
        exit()
fh_cand_SDP.close()

#Check that each candidate SDP directory was successfully created and contains fasta-files for processing

SDP_dirs = cand_SDP_dict.keys()
SDP_dirs_verified = list()
run_dir = os.getcwd()
for dir in SDP_dirs:
    isdir = os.path.isdir(dir)
    if (isdir == True):
        os.chdir(dir)
        ffn_suffix = '*ffn'
        ffn_files = glob.glob(ffn_suffix)
        try:
            open(ffn_files[0])
            SDP_dirs_verified.append(dir)
        except Exception:
            pass
        os.chdir(run_dir)

#Print verified SDP-directories as string to stdout for further processing with bash script

if (len(SDP_dirs_verified) > 0):
    SDP_dirs_verified.sort()
    SDP_dir_string = " ".join(SDP_dirs_verified)
    print(SDP_dir_string)
