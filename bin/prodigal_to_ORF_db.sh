#!/bin/bash

#Parse script arguments (name of ortholog-file and output directory), or print script usage if required arguments were not provided
usage="
bash $(basename "$0") [-h] [-c THRESHOLD] -f [file_list] [-o OUT_DB]

Generates a database of metagenomic ORFs from ORFs predicted with prodigal on multiple independent metagenomic assemblies. Sample-names will be added to the fasta-headers, based on the input fasta file-names, so that each header is unique in the database. Clustering to remove redundancy across assemblies can optionally be done, with cd-hit. If clustering, a threshold of 99% is recommended (-c 99)

The following options can/should be provided:
    -f file-list; simple text-file with the names of files to be included in the db (required)
    -o output; name for the output database
    -c threshold of percentage identity used for clustering with cd-hit (optional)"

options=':hc:f:o:'
while getopts $options option; do
    case "${option}" in
        h) echo "$usage"; exit;;
        f) FILE_LIST=${OPTARG};;
        o) OUT_DB=${OPTARG};;
        c) THRESHOLD=${OPTARG};;
       \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
    esac
done

if [ ! "$FILE_LIST" ] ; then
    echo "ERROR: Please provide a text-file listing the input-files for the database"
    echo "$usage" >&2; exit 1
fi
if [ ! "$OUT_DB" ] ; then
    echo "ERROR: Please specify a name for the output database"
    echo "$usage" >&2; exit 1
fi

#Read the input file-list, and edit the headers for each file

rm $OUT_DB
while IFS= read -r line; do
    echo "Processing input-file: $line"
    python3 edit_ORF_headers.py $line >> $OUT_DB
done < "$FILE_LIST"
echo "ORF-db written to file: $OUT_DB"