#!/usr/bin/env bash

#Parse script arguments (name of ortholog-file and output directory), or print script usage if required arguments were not provided

usage="
bash $(basename "$0") [-h] [-c CAND_SDP] [-i IN_DIR] [-d ORF_DB]

Perform a validation analysis of candidate SDPs (inferred from sequenced isolates), with metagenomic data. Candidate SDPs are detailed in the input SDP file (-c CAND_SDP), with genome-id in column1 and SDP affiliation in column three. 

The script will perform the following four steps:

1. Core gene sequences of the candidate SDPs are subset from fasta-files in the user-provided input directory (-i IN_DIR), and written into new SDP directories
2. For each SDP, the core gene sequences are blasted against the metagenomic database (-d ORF_DB)
3. Significant blast-hits to ORFs (min-length 200bp) are printed to the SDP directories
4. The recruited ORFs are added to the alignments in the input directory (IN_DIR), one at a time, and the maximum percentage identity within the alignment is recorded (and the SDP affiliation)

The following options should be provided:

-c: CAND_SDP. Simple tab-delimited file with genome-id in tab1 and SDP affiliation in tab3 (required)
-i: IN_DIR. Input-directory, containing fasta-files (*ffn) and the corresponding untrimmed alignments (*aln_nuc.fasta) (required)
-d: ORF_DB. Name of ORF database file from which ORFs should be recruited (or full path, if the file is not in the run-directory) (required)
"

options=':hc:i:d:'
while getopts $options option; do
    case "${option}" in
        h) echo "$usage"; exit;;
        c) CAND_SDP=${OPTARG};;
        i) IN_DIR=${OPTARG};;
        d) ORF_DB=${OPTARG};;
       \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
    esac
done

if [ ! "CAND_SDP" ] ; then
    echo "ERROR: Please provide a tab-delimited file detailing genome-ids for candidate SDPs"
    echo "$usage" >&2; exit 1
fi
if [ ! "$IN_DIR" ] ; then
    echo "ERROR: Please specify the input directory, containing the alignment files that should be subset"
    echo "$usage" >&2; exit 1
fi
if [ ! "$ORF_DB" ] ; then
    echo "ERROR: Please provide the name of the database, from which ORFs should be recruited"
    echo "$usage" >&2; exit 1
fi

#STEP1: Subset SDP-data from fasta-files in the input directory. Generate SDP directories, and write subset data to the corresponding dirs.

echo 'Performing step1: subset SDP-data fasta-files in the input directory' | tee -a "log.txt"
python3 bin/extract_SDP_coreseqs.py $CAND_SDP $IN_DIR

#STEP2: Blasting SDP core sequences against ORF database

SDP_DIR_STRING=$(python3 bin/get_SDP_dirs.py $CAND_SDP) #Get the SDP-dirs to be processed from the candidate-SDP input file, check that they contain fasta-files for blasting (exit bash-script with error message if they dont)
IFS=' ' read -ra SDP_DIR <<< $SDP_DIR_STRING
if [ ${SDP_DIR[0]} == "ERROR" ]; then
    echo $SDP_DIR_STRING | tee -a "log.txt"
    exit
fi
echo 'Performing step2: blasting SDP core sequences against ORF database' | tee -a "log.txt"
RUN_DIR=$(pwd)
for DIR in "${SDP_DIR[@]}"; do
    echo "Blasting core sequences in directory: "$DIR 
    cd $DIR 
    COUNTER=0
    for i in $(ls *ffn); do 
        BLAST_OUTFILE=${i:0:9}".blastn"
        if [[ -f "$BLAST_OUTFILE" ]]; then
            echo "This blast-file already exists: "$BLAST_OUTFILE" in "$DIR | tee -a "log.txt"
        else
            blastn -db $RUN_DIR/$ORF_DB -query $i -outfmt 5 -evalue 1e-5 -perc_identity 70 > $BLAST_OUTFILE
        fi
        (( COUNTER++ ))       
        if (( $COUNTER % 10 == 0 ))
        then
            echo "Finished blasting $COUNTER files.."
        fi
    done    
    cd $RUN_DIR
    echo "Done" 
done
    
#STEP3: Recruiting ORF sequences from ORF-db, based on blast-files

echo 'Performing step3: recruiting ORF sequences from ORF database file, based on blast-files' | tee -a "log.txt"
for DIR in "${SDP_DIR[@]}"; do
    echo "Processing SDP-directory: "$DIR 
    python3 bin/recruit_ORFs.py $ORF_DB $DIR
done

#STEP4: Adding recruited ORFs to core alignments, recording max perc-id and SDP affiliation. 

echo 'Performing step4: adding recruited ORFs to core sequence gene alignments, and calculating max percentage identity' | tee -a "log.txt"
for DIR in "${SDP_DIR[@]}"; do
    echo "Processing SDP-directory: "$DIR | tee -a "log.txt"
    python3 bin/orf_aln_perc_id.py $CAND_SDP $IN_DIR $DIR
done 

echo 'All done! Check log.txt file for some summary stats on the results'
























