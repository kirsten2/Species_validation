#!/usr/bin/env python3
import sys

#Usage: python3 get_scp_orthologs.py OrthologousGroups.txt 

try:
    ortho_in = sys.argv[1]
except:
    print('Provide a text-file with orthologous gene-families, formatted as per orthofinder. Exiting script')
    exit()

def genome_count_line(line):
    genomes = dict()
    split_line = line.split()
    OG_id = split_line.pop(0)
    for seq_id in split_line:
        split_id = seq_id.split('_')
        genome_id = split_id[0]
        genomes[genome_id] = genomes.get(genome_id,0) + 1
    return(genomes)
        
        
#Get all the genome-identifiers present in the ortholog-file
genome_ids = dict()
fh_ortho_in = open(ortho_in)
for line in fh_ortho_in:
    line = line.strip()
    genome_count = genome_count_line(line)
    for genome in genome_count:
        genome_ids[genome] = 1
fh_ortho_in.close()
nb_genome_ids = len(genome_ids)

#Go through ortholog-file again, and check whether all genome-ids are present in a single copy, for each family. If so, print to file.

fh_ortho_in = open(ortho_in)
fh_out = open('single_ortho.txt','w')
for line in fh_ortho_in:
    line = line.strip()
    genome_count = genome_count_line(line)
    all_counts_equals_one = all(x == 1 for x in genome_count.values())
    if (len(genome_count) == nb_genome_ids and all_counts_equals_one == True):
        fh_out.write(line + '\n')
fh_ortho_in.close()
fh_out.close()

  
    
