#!/usr/bin/env python3
import sys
import os
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SearchIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from shutil import copyfile
import statistics

#Usage: python3 orf_aln_perc_id.py Candidate_species.txt core_aln_dir_name species_dir

def get_cand_species_dict(cand_species_file):
    fh_cand_species_file = open(cand_species_file)
    cand_species_dict = dict()
    for line in fh_cand_species_file:
        line = line.strip()
        if (line.find('Genome_id') != -1): continue
        split_line = line.split('\t')
        genome_id = split_line[0]
        species = split_line[2]
        cand_species_dict[genome_id] = species
    return(cand_species_dict)  

def perc_id(aln, orf_length):
    nb_col = aln.get_alignment_length()
    nb_col_ungapped = 0
    nb_matches = 0
    i = 0
    perc_id = 0
    while(i < nb_col): #iterate over all alignment columns
        aln_col = aln[:,i]
        if (aln_col.find('-') == -1): #no gap for this column
            nb_col_ungapped += 1 
            if (aln_col[0] == aln_col[1]):
                nb_matches += 1
        i += 1
    if (nb_col_ungapped/orf_length > 0.5): #only calculate perc-id if at least half of the ORF aligns to the core gene sequence, otherwise return perc_id value of zero
        perc_id = (nb_matches/nb_col_ungapped)*100
    return(perc_id)

def add_orf_perc_id(ref_aln_file, orf_file):
    #Add orf to reference alignment with muscle
    process = subprocess.Popen(['muscle', '-profile', '-in1',ref_aln_file, '-in2', orf_file],
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = process.communicate()
    new_aln =  AlignIO.read(StringIO(stdout), "fasta")
    #Get seq-ids and corresponding seq-objects for the reference genomes and orf
    ref_aln_obj = AlignIO.read(ref_aln_file,"fasta")
    ref_aln_ids = [aln_record.id for aln_record in ref_aln_obj]
    all_seq_objects = dict()
    for aln_record in new_aln:
        all_seq_objects[aln_record.id]=aln_record
    orf_seq_obj = SeqIO.read(orf_file, "fasta")
    orf_length = len(orf_seq_obj.seq)
    orf_id = orf_seq_obj.id
    #Generate pairwise alignments by subsetting the new core alignment containing the orf. For each pairwise alignment, calculate the percentage identity. Save the maximum percentage id and corresponding reference genome id
    max_perc_id = 0
    max_perc_id_gene = None
    for ref_id in ref_aln_ids:
        ref_id_split = ref_id.split('_')
        if (ref_id_split[0] in species_dict and species_dict[ref_id_split[0]] == species_dir): #Only align orf to core-seqs of the recruiting species
            sub_seq_objects = [all_seq_objects[ref_id], all_seq_objects[orf_id]]
            sub_aln = MultipleSeqAlignment(sub_seq_objects)
            sub_aln_perc_id = perc_id(sub_aln, orf_length) 
            if (sub_aln_perc_id > max_perc_id and sub_aln_perc_id > 60): 
                max_perc_id = sub_aln_perc_id
                max_perc_id_gene = ref_id
    return(orf_id, max_perc_id_gene, max_perc_id)
 
def get_best_blast_hit(orf_file, db_file):
    blastn_cmd = NcbiblastnCommandline(query=orf_file, db=db_file, evalue=0.001, outfmt=5)
    stdout, stderr = blastn_cmd()
    blast_obj =  SearchIO.read(StringIO(stdout), "blast-xml")
    top_hit_genome = None
    if (len(blast_obj) != 0):
        top_hit_id = blast_obj[0].id
        split_hit_id = top_hit_id.split('_')
        top_hit_genome = split_hit_id[0]
    return(top_hit_genome)
      
#Parse input arguments and working dirs
fh_log_out = open('log.txt','a')
species_dict = get_cand_species_dict(sys.argv[1])
current_dir = os.getcwd()
core_aln_dir = current_dir + '/' + sys.argv[2]
species_dir = sys.argv[3]

#Change into species-dir, use glob to list ORF-files to be processed (remove previous perc-id file if present)
os.chdir(species_dir)
if os.path.isfile('perc_id.txt'):
    os.remove('perc_id.txt')
orf_file_suffix = '*_orfs.ffn'
orf_files = glob.glob(orf_file_suffix)
orf_files.sort()
nb_orf_files = len(orf_files)
if (nb_orf_files == 0):
    os.chdir(current.dir)
    fh_log_out.write('ERROR: No orf-files with the suffix "*_orfs.ffn" were found in the species directory. Check the directory: \n')
    fh_log_out.write(species_dir + "\n")
    fh_log_out.close()
    exit()

#Loop over the orf-files in the species dir, processing them one by one
fh_perc_id_out = open('perc_id.txt','a')
header = "\t".join(["OG_group", "ORF_id","Gene_id","Perc_id","Closest_SDP"])
fh_perc_id_out.write(header + '\n')
count_aln_results = 0
count_species_recruits = dict()
for orf_file in orf_files:
    print('Processing orf-file: ', orf_file)
    filename_split = orf_file.split('_')
    OG = filename_split[0]
    core_aln_file = OG + '_aln_nuc.fasta'
    core_aln_file_fullname = core_aln_dir + '/' + core_aln_file
    core_ffn_file = OG + '.ffn'
    core_ffn_file_full = core_aln_dir + '/' + core_ffn_file
    #copy core ffn-file as temporary file to SDP-dir, makeblastdb
    copyfile(core_ffn_file_full, "temp.ffn")
    makeblastdb_cmd = NcbimakeblastdbCommandline(dbtype="nucl", input_file="temp.ffn")
    makeblastdb_cmd()
    #Loop over orfs in current orf-file. Write each orf to temporary fasta-file, blast against the core-seqs and get species-affiliation for first blast-hit. 
    for seq_record in SeqIO.parse(orf_file, "fasta"):
        SeqIO.write(seq_record, "temp_orf.ffn", "fasta")
        blast_result = get_best_blast_hit('temp_orf.ffn', 'temp.ffn')
        if (blast_result == None): continue
        #Add  orf to core alignment with muscle and get max perc-id
        aln_result = add_orf_perc_id(core_aln_file_fullname, 'temp_orf.ffn')       
        if (aln_result[2] != 0):
            orf_max_perc_id = round(aln_result[2],2)
            orf_id = aln_result[0]
            best_hit_to_species_gene_id = aln_result[1]
            best_species = "other"
            if blast_result in species_dict:
                best_species = species_dict[blast_result]
            if (best_species == species_dir):
                count_species_recruits[orf_id] = orf_max_perc_id
            out_str = "\t".join([OG, orf_id, best_hit_to_species_gene_id, str(orf_max_perc_id), best_species])
            fh_perc_id_out.write(out_str + "\n")
            count_aln_results += 1
fh_perc_id_out.close

#Final results and clean-up
temp_files = glob.glob('temp*')
for file in temp_files:
    os.remove(file)
os.chdir(current_dir)
fh_log_out.write('A total of ' + str(count_aln_results) + ' ORFs were aligned, and printed to file "perc_id.txt"\n')
nb_species_recruits = len(count_species_recruits)
if (nb_species_recruits == 0):
    fh_log_out.write('However, none of these had a best blast hit to the recruiting species\n\n')
else:
    species_recruit_values = list(count_species_recruits.values())
    median_species_perc_id = species_recruit_values[0]
    if (nb_species_recruits > 1):
        median_species_perc_id = statistics.median(species_recruit_values)
    fh_log_out.write('Out of these, ' + str(nb_species_recruits) + ' had a best blast hit to the recruiting species, with a median perc-id of ' + str(round(median_species_perc_id,2)) + "\n\n")
fh_log_out.close()
