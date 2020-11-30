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
import numpy as np
from shutil import copyfile
import statistics

#Usage: python3 orf_aln_perc_id.py Candidate_SDPs.txt core_aln_dir_name SDP_dir

def get_cand_SDP_dict(cand_SDP_file):
    fh_cand_SDP_file = open(cand_SDP_file)
    cand_SDP_dict = dict()
    for line in fh_cand_SDP_file:
        line = line.strip()
        if (line.find('Genome_id') != -1): continue
        split_line = line.split('\t')
        genome_id = split_line[0]
        SDP = split_line[2]
        cand_SDP_dict[genome_id] = SDP
    return(cand_SDP_dict)  
    
def trim_aln(aln):
    nb_seqs = len(aln)
    nb_col = aln.get_alignment_length()
    seq_ids = [aln_record.id for aln_record in aln]
    #Go through alignment columns, and select which ones to keep (gap in less than 50% of the sequences)
    i = 0
    gap_char = '-'
    selected_columns = list()
    while (i < nb_col):
        aln_col = aln[:,i]
        nb_gaps = aln_col.count(gap_char)
        gap_fraction = round((nb_gaps/nb_seqs),3)
        if (gap_fraction < 0.5):
            selected_columns.append(i)
        i += 1
    #Convert alignment to numpy object, create trimmed alignment as new numpy object
    aln_np = np.array([list(rec) for rec in aln], dtype=str)
    nb_seqs = len(aln)
    aln_np_trim = np.empty((nb_seqs,0),dtype=str)
    for j in selected_columns:
        col_seq = aln_np[:,j].tolist()
        aln_np_trim = np.append(aln_np_trim,np.array([col_seq]).transpose(),axis=1)
    #Convert trimmed numpy object to alignment object
    k = 0
    seq_records = list()
    while(k < nb_seqs):
        row_seq = aln_np_trim[k,:].tolist()
        row_seq_str = ''.join(row_seq)
        row_seq_id = seq_ids[k]
        new_seq_record = SeqRecord(Seq(row_seq_str), id = row_seq_id, description="")
        seq_records.append(new_seq_record)
        k += 1
    new_aln = MultipleSeqAlignment(seq_records)
    return(new_aln)

def perc_id(aln, orf_length):
    nb_col = aln.get_alignment_length()
    nb_col_ungapped = 0
    nb_matches = 0
    i = 0
    perc_id = 0
    while(i < nb_col):
        aln_col = aln[:,i]
        if (aln_col.find('-') == -1): #no gap..
            nb_col_ungapped += 1 
            if (aln_col[0] == aln_col[1]):
                nb_matches += 1
        i += 1
    if (nb_col_ungapped/orf_length > 0.5):
        perc_id = (nb_matches/nb_col_ungapped)*100
    return(perc_id)

def add_orf_perc_id(ref_aln_file, orf_file):
    #Add orf to reference alignment with muscle, trim the alignment
    process = subprocess.Popen(['muscle', '-profile', '-in1',ref_aln_file, '-in2', orf_file],
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = process.communicate()
    new_aln =  AlignIO.read(StringIO(stdout), "fasta")
    trimmed_aln = trim_aln(new_aln)
    #Get seq-ids and corresponding seq-objects for the reference genomes and orf
    ref_aln_obj = AlignIO.read(ref_aln_file,"fasta")
    ref_aln_ids = [aln_record.id for aln_record in ref_aln_obj]
    all_seq_objects = dict()
    for aln_record in trimmed_aln:
        all_seq_objects[aln_record.id]=aln_record
    orf_seq_obj = SeqIO.read(orf_file, "fasta")
    orf_length = len(orf_seq_obj.seq)
    orf_id = orf_seq_obj.id
    #Generate pairwise alignments by subsetting the trimmed alignment. For each pairwise alignment, calculate the percentage identity. Save the maximum percentage id and corresponding reference genome id
    max_perc_id = 0
    max_perc_id_genome = None
    for ref_id in ref_aln_ids:
        ref_id_split = ref_id.split('_')
        if (SDP_dict[ref_id_split[0]] == SDP_dir): #Only align to core-seqs of recruiting SDP
            sub_seq_objects = [all_seq_objects[ref_id], all_seq_objects[orf_id]]
            sub_aln = MultipleSeqAlignment(sub_seq_objects)
            sub_aln_perc_id = perc_id(sub_aln, orf_length) 
            if (sub_aln_perc_id > max_perc_id): 
                max_perc_id = sub_aln_perc_id
                max_perc_id_genome = ref_id
    return(orf_id, max_perc_id_genome, max_perc_id)
 
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
SDP_dict = get_cand_SDP_dict(sys.argv[1])
current_dir = os.getcwd()
core_aln_dir = current_dir + '/' + sys.argv[2]
SDP_dir = sys.argv[3]

#Change into SDP-dir, use glob to list ORF-files to be processed (remove previous perc-id file if present)
os.chdir(SDP_dir)
if os.path.isfile('perc_id.txt'):
    os.remove('perc_id.txt')
orf_file_suffix = '*_orfs.ffn'
orf_files = glob.glob(orf_file_suffix)
orf_files.sort()
nb_orf_files = len(orf_files)
if (nb_orf_files == 0):
    os.chdir(current.dir)
    fh_log_out.write('ERROR: No orf-files with the suffix "*_orfs.ffn" were found in the SDP directory. Check the directory: \n')
    fh_log_out.write(sdp_dir + "\n")
    fh_log_out.close()
    exit()

#Loop over the orf-files, get percentage id and print to file
fh_perc_id_out = open('perc_id.txt','a')
header = "\t".join(["OG_group", "ORF_id","Gene_id","Perc_id","Closest_SDP"])
fh_perc_id_out.write(header + '\n')
count_aln_results = 0
count_SDP_recruits = dict()
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
    #Write each orf to temporary fasta-file, blast against core-seqs, add to core alignment with muscle and get max perc-id
    for seq_record in SeqIO.parse(orf_file, "fasta"):
        SeqIO.write(seq_record, "temp_orf.ffn", "fasta")
        blast_result = get_best_blast_hit('temp_orf.ffn', 'temp.ffn')
        if (blast_result == None): continue
        aln_result = add_orf_perc_id(core_aln_file_fullname, 'temp_orf.ffn')       
        if (aln_result[2] != 0):
            orf_max_perc_id = round(aln_result[2],2)
            orf_id = aln_result[0]
            best_hit_to_SDP_gene_id = aln_result[1]
            best_SDP = SDP_dict[blast_result]
            if (best_SDP == SDP_dir):
                count_SDP_recruits[orf_id] = orf_max_perc_id
            out_str = "\t".join([OG, orf_id, best_hit_to_SDP_gene_id, str(orf_max_perc_id), best_SDP])
            fh_perc_id_out.write(out_str + "\n")
            count_aln_results += 1
fh_perc_id_out.close

#Final results and clean-up
temp_files = glob.glob('temp*')
for file in temp_files:
    os.remove(file)
os.chdir(current_dir)
fh_log_out.write('A total of ' + str(count_aln_results) + ' ORFs were aligned, and printed to file "perc_id.txt"\n')
nb_SDP_recruits = len(count_SDP_recruits)
SDP_recruit_values = list(count_SDP_recruits.values())
mean_SDP_perc_id = statistics.mean(SDP_recruit_values)
mean_SDP_perc_id_round = round(mean_SDP_perc_id, 2) 
fh_log_out.write('Out of these,' + str(nb_SDP_recruits) + ' had a best blast hit to the recruiting SDP, with a mean perc-id of ' + str(mean_SDP_perc_id_round) + "\n\n")
fh_log_out.close()
