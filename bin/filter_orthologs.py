#!/usr/bin/env python3
import sys
import os
import os.path
from Bio import SeqIO
from Bio.Seq import Seq

#Usage: python3 filter_orthologs.py single_ortho.txt

def get_OG_fams(file):
    OG_fams = dict()
    fh_file = open(file)
    for line in fh_file:
        line = line.strip()
        split_line = line.split()
        OG_id = split_line.pop(0)[:-1]
        OG_fams[OG_id] = split_line
    return(OG_fams)

def get_genome_ids(seq_ids):
    genomes = list()
    for seq_id in seq_ids:
        split_id = seq_id.split('_')
        genome_id = split_id[0]
        genomes.append(genome_id)
    return(genomes)
    
def get_gene_lengths(genome_id, genome_dir):
    ffn_file = genome_dir + genome_id + '.ffn'
    seq_length_dict = dict()
    for seq_record in SeqIO.parse(ffn_file, "fasta"):
        seq_length = len(seq_record.seq)
        seq_length_dict[seq_record.id] = seq_length
    return(seq_length_dict)

def count_small_seqs(OG_id):
    small_seq_count = 0
    OG_gene_ids = OG_fams[OG_id]
    for gene in OG_gene_ids:
        if (all_gene_lengths[gene] < 300):
            small_seq_count += 1
    return(small_seq_count)
            

#Check that the input-file "single_ortho.txt" has been provided. Determine full path to ffn-files
try:
    ortho_in = sys.argv[1]
    fh_ortho_in = open(ortho_in)
except:
    print('The file "single_ortho.txt", was not found in the run-directory. Exiting script')
    exit()
current_dir = os.getcwd()
genome_ffn_dir = current_dir + '/genome_ffn/'

#Get all the genome-ids from the ortholog-file, and all the gene-ids associated with each gene-family
genome_ids = dict()
OG_fams = dict()
for line in fh_ortho_in:
    line = line.strip()
    split_line = line.split()
    OG_id = split_line.pop(0)[:-1]
    OG_fams[OG_id] = split_line
    genome_ids_OG = get_genome_ids(split_line)
    for genome in genome_ids_OG:
        genome_ids[genome] = 1
fh_ortho_in.close()

#Get the gene-lengths of all genes, for each genome
all_gene_lengths = dict()
for genome in genome_ids.keys():
    gene_lengths = get_gene_lengths(genome, genome_ffn_dir)
    all_gene_lengths.update(gene_lengths)

#Generate the filtered ortholog-file
fh_out = open('single_ortho_filt.txt','w')
OG_fam_ids = list(OG_fams.keys())
OG_fam_ids.sort()
for id in OG_fam_ids:
    small_seq_count_OG = count_small_seqs(id)
    if (small_seq_count_OG == 0):
        OG_fam_str = ' '.join(OG_fams[id])
        str_out = id + ': ' + OG_fam_str
        fh_out.write(str_out + '\n')
fh_out.close()
