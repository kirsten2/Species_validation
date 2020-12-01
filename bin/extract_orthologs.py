#!/usr/bin/env python3
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#Usage: python3 extract_orthologs.py firm5

def genome_count_line(seq_ids):
    genomes = dict()
    for seq_id in seq_ids:
        split_id = seq_id.split('_')
        genome_id = split_id[0]
        genomes[genome_id] = genomes.get(genome_id,0) + 1
    return(genomes)

def get_seq_objects(filename):
    seq_record_dict = dict()
    for seq_record in SeqIO.parse(filename, "fasta"):
        if (seq_record.seq[-1] == '*'):  #Sometimes, stop-codons are indicated with an asterisk for amino acid gene sequences, which in turn generates warnings when using "muscle" for alignments. This block of code remove the trailing asterisk when present.
            seq_string = str(seq_record.seq)
            seq_chomp = seq_string[:-1]
            new_record =  SeqRecord(
                Seq(seq_chomp),
                id=seq_record.id,
                name="",
                description=""
                )
            seq_record_dict[seq_record.id] = new_record
        else:
            seq_record_dict[seq_record.id] = seq_record
    return(seq_record_dict)

#Open the ortholog file
try:
    fh_ortho_in = open('single_ortho_filt.txt')
except:
    print('Input-file not found: "single_ortho_filt.txt". Exiting script')
    exit()

#Get all the genome-ids present in the ortholog-file, and all the gene-ids associated with each gene-family
genome_ids = dict()
OG_fams = dict()
for line in fh_ortho_in:
    line = line.strip()
    split_line = line.split()
    OG_id = split_line.pop(0)[:-1]
    OG_fams[OG_id] = split_line
    genome_count = genome_count_line(split_line)
    for genome in genome_count:
        genome_ids[genome] = 1
fh_ortho_in.close()

#Construct dictionaries of gene-sequences 
current_dir = os.getcwd()
genome_faa_dir = current_dir + '/genome_faa/'
genome_ffn_dir = current_dir + '/genome_ffn/'
ffn_seq_objects = dict()
faa_seq_objects = dict()
for genome in genome_ids.keys():
    ffn_file = genome_ffn_dir + genome + '.ffn'
    faa_file = genome_faa_dir + genome + '.faa'
    try:
        fh_ffn_file = open(ffn_file)
        fh_faa_file = open(faa_file)
    except:
        print('One or both of these files are missing: ')
        print(ffn_file)
        print(faa_file)
        print('Exiting script')
        exit()
    genome_ffn_seq_objects = get_seq_objects(ffn_file)
    ffn_seq_objects.update(genome_ffn_seq_objects)
    genome_faa_seq_objects = get_seq_objects(faa_file)
    faa_seq_objects.update(genome_faa_seq_objects)


#Create the output directory, print multi-fasta files corresponding to each gene-family into the directory
output_dir_name = sys.argv[1]
output_dir = current_dir + '/' + output_dir_name
try:
    os.mkdir(output_dir)
except OSError:
    print('Could not create the output dir:')
    print(output_dir)
    exit()
os.chdir(output_dir)

#Print multi-fasta files corresponding to each gene-family in the output dir
for OG in OG_fams.keys():
    OG_seq_ids = OG_fams[OG]
    OG_ffn_seq_obj = [ffn_seq_objects[x] for x in OG_seq_ids]
    OG_faa_seq_obj = [faa_seq_objects[x] for x in OG_seq_ids]
    ffn_outfile = OG + '.ffn'
    faa_outfile = OG + '.faa'
    SeqIO.write(OG_ffn_seq_obj,ffn_outfile,"fasta")
    SeqIO.write(OG_faa_seq_obj,faa_outfile,"fasta")
    

    
