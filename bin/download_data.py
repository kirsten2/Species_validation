#!/usr/bin/env python

# ============================================================================ #
# download_db.py: download example data for species-validation after cloning from github
#
# Author: Kirsten Ellegaard (kirsten.ellegaard@unil.ch)
#
# ============================================================================ #

import os
import sys
import tempfile
import shutil
import subprocess
import hashlib
import time
import urllib.request
import argparse

#Links, file-names, sizes and md5 check-sums
zenodo_links = {
    "genes": "https://sandbox.zenodo.org/record/709842/files/genes_db.tar.gz",
    "genomes": "https://sandbox.zenodo.org/record/709842/files/genomes_db.tar.gz",
    "metagenomes": "https://sandbox.zenodo.org/record/709842/files/metagenomic_orfs.tar.gz"
}
file_names = {
    "genes": "genes_db.tar.gz",
    "genomes": "genomes_db.tar.gz",
    "metagenomes": "metagenomic_orfs.tar.gz"
}
file_sizes = {
    "genes": 130,
    "genomes": 96,
    "metagenomes": 43
}
check_sums = {
    "genes": "fee84a9017e8dc8b1630991cb05ca787",
    "genomes": "ed010b67bdcb1fa931b7de16aea169d8", 
    "metagenomes": "21a985aa90bc49f866c65883431116c7"
}

#function to print progress bar for python 3
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r %d%%, %d MB, %d KB/s, %d seconds passed" % (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

#function to download data
def save_f(url, filename):
    if "--no-download-progress" in sys.argv:
        urllib.request.urlretrieve(url, filename)
    else:
        urllib.request.urlretrieve(url, filename, reporthook)

#function to calculate md5
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

#Parse input arguments, check for existing files/directories

parser = argparse.ArgumentParser(description="This script downloads data from zenodo, for running various metagenomic pipelines. Specify at least one data-set from the following options:")
parser.add_argument("--genes", help="Download the gene database", action="store_true")
parser.add_argument("--genomes", help="Download the genome database", action="store_true")
parser.add_argument("--metagenomes", help="Download metagenomic ORF example data", action="store_true")
args = parser.parse_args()

current_dir = os.getcwd()
if args.genes:
    genes_faa_dir = current_dir + '/genes_faa/'
    genes_ffn_dir = current_dir + '/genes_ffn/'
    if (os.path.isdir(genes_faa_dir) or os.path.isdir(genes_ffn_dir)):
        print("One or both of these directories already exist:")
        print(genes_faa_dir)
        print(genes_ffn_dir)
        print("Exiting script!")
        sys.exit(1)
if args.genomes:
    bed_dir = current_dir + '/bed_files/'
    db_file = 'genomes_db.fasta'
    if (os.path.isdir(bed_dir)):
        print("Directory with bed-files already exist:")
        print(bed_dir)
        print("Exiting script!")
        sys.exit(1)
    if (os.path.isfile(db_file)):
        print('The file "genomes_db.fasta" already exist in the run directory, exiting script!')
        sys.exit(1)
if args.metagenomes:
    meta_orf_file = 'metagenomic_orfs.ffn'
    if (os.path.isfile(meta_orf_file)):
        print('The file "metagenomic_orfs.ffn" already exist in the run directory, exiting script!')
        sys.exit(1)
if not any(vars(args).values()):
    parser.print_help(sys.stderr)
    sys.exit(1)

#Download compressed files, check md5 and unpack

all_args = vars(args)
for arg in all_args.keys():
    if (all_args[arg] == False): continue
    #Download compressed file
    print('Downloading gzipped file:', file_names[arg], '(~' + str(file_sizes[arg]) +'Mb)')
    save_f(zenodo_links[arg], file_names[arg])
    #Checking md5
    print('\nChecking MD5..')
    download_md5 = md5(file_names[arg])
    if (download_md5 == check_sums[arg]):
        print('OK!')
    else:
        print('MD5 verification failed for', file_names[arg])
        print('Removing corrupt file, and exiting script')
        os.remove(file_names[arg])
    #Extract files from archive
    print('Extracting files from:', file_names[arg])
    extract_cmd = "tar -zxvf " + file_names[arg] + " -C "+ current_dir
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        print("Error: failed to extract files\n")
        exit()
    if process.returncode:
        print("Error: failed to extract files\n")
        exit()
    os.remove(file_names[arg])
    
print('All done!')

