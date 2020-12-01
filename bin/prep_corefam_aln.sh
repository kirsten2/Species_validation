#!/bin/bash

##Parse script arguments (name of ortholog-file and output directory), or print script usage if they were not provided
usage="
bash bin/$(basename "$0") [-h] [-d OUTDIR] [-o ORTHOFILE]

Generates back-translated nucleotide alignments of filtered single-copy orthologous gene-families, to be used for phylogenetic inference and/or SDP metagenomic validation pipeline. 

The following options must be provided:
     -d name of directory for alignment files (required)
     -o name of orthofinder file, from which single-copy families will be parsed (required)"

options=':ho:d:'
while getopts $options option; do
    case "${option}" in
        h) echo "$usage"; exit;;
        o) ORTHOFILE=${OPTARG};;
        d) OUTDIR=${OPTARG};;
        :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
       \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
    esac
done
if [ ! "$OUTDIR" ] || [ ! "$ORTHOFILE" ]; then
    echo "$usage" >&2; exit 1
fi

##Check script requirements

RUN_DIR=$(pwd)

#Check that the output directory name doesnt already exist
if [ -d $OUTDIR ]; then
    echo "Output directory already exists! Exiting script"
    exit 1
fi

#Check that input orthofile is present in the run directory
if [[ ! -f "$ORTHOFILE" ]]; then
    echo "The orthofinder file could not be opened:"
    echo $ORTHOFILE
    echo "Exiting script"
    exit 1
fi

#Check that muscle alignment software is installed and in path
if ! [ -x "$(command -v muscle)" ]; then
  echo 'Error: muscle is not installed/not in path.' >&2
  exit 1
fi

#Check that the directories "genome_faa" and "genome_ffn" are present in the run-directory

if [ ! -d "genome_faa" ] || [ ! -d "genome_ffn" ]; then
    echo 'the directory "genome_faa"  and/or "genome_ffn" doesnt exist in the run directory:'
    echo $RUN_DIR
    echo "Exiting script"
    exit 1
fi

##Start the script

#Generate ortholog-file with single-copy ortholog families. Filter single-copy orthologs by length (gene-families containing genes shorter than 300bp are removed)

python3 bin/get_scp_orthologs.py $ORTHOFILE 
python3 bin/filter_orthologs.py single_ortho.txt

#Extract gene-sequences for each single-copy ortholog family into fasta-files in the output directory (exit script if unsuccessful).

echo "Extracting filtered single-copy ortholog gene-family sequences into directory: $OUTDIR"
python3 bin/extract_orthologs.py $OUTDIR
if [ ! -d $OUTDIR ]; then
    exit 1 
fi

#For each gene-family: 1. align the amino-acid sequences, 2. translate to codon-aligned nucleotide alignement, 3. simplify seq-headers to genome-id

cd $OUTDIR
echo "Starting alignment of amino-acid sequences, and back-translation to nucleotide alignments"
SAVEIFS=$IFS  #code to counter-act automatic backups on mac, generating file-names with spaces..
IFS=$(echo -en "\n\b")
for i in $( ls *.faa ); do
    if [[ ! "$i" =~ [[:space:]] ]]; then 
        OG=${i:0:9}
        echo "Processing: "$OG
        muscle -in $OG".faa" -out $OG"_aln.fasta" -quiet
        #mafft --auto --quiet $OG".faa" > $OG"_aln.fasta"
        ALN_AA=$OG"_aln.fasta"
        FFN_FILE=$OG".ffn"
        python3 ../bin/aln_aa_to_dna.py $ALN_AA $FFN_FILE
        ALN_NUC=$OG"_aln_nuc.fasta"
        sed 's/_.*$//g' $ALN_NUC > temp_aln.txt
        mv temp_aln.txt $ALN_NUC
    fi
done
IFS=$SAVEIFS
echo "All done!"