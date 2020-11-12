#!/bin/bash

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

echo 'Performing step1: subset SDP-data fasta-files in the input directory'
python3 bin/extract_SDP_coreseqs.py $CAND_SDP $IN_DIR

#STEP2: Blasting SDP core sequences against ORF database
RUN_DIR=$(pwd)
echo 'Performing step2: blasting SDP core sequences against ORF database'
while IFS= read -r DIR; do
    cd $DIR
    echo "Blasting core sequences in directory:"
    echo $DIR
    COUNTER=0
    for i in $(ls *ffn); do 
        OG=${i:0:9}
        blastn -db $RUN_DIR/$ORF_DB -query $i -outfmt 5 -evalue 1e-5 -perc_identity 70 > $OG".blastn"
        (( COUNTER++ ))       
        if (( $COUNTER % 10 == 0 ))
        then
            echo "Finished blasting $COUNTER files.."
        fi
    done    
    cd $RUN_DIR
done < "SDP_dirs.txt"


#STEP3: Extracting ORF sequences from ORF-db, based on blast-files, print to SDP-dirs

cd $RUN_DIR
while IFS= read -r DIR; do
    echo "Parsing blast-files from directory $DIR"
    python3 bin/recruit_ORFs.py $ORF_DB $DIR
done < "SDP_dirs.txt"

#STEP4: TO BE DONE. Add recruited ORFs to core alignments, record max perc-id and SDP affiliation. Print to file in SDP_dir




































