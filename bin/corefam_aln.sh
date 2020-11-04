#!/bin/bash

#Parse script arguments (name of ortholog-file and output directory), or print script usage if they were not provided
usage="
bash bin/$(basename "$0") [-h] [-d OUTDIR] [-o ORTHOFILE]

Generates trimmed nucleotide alignments of filtered single-copy orthologous gene-families, to be used for phylogenetic inference and/or SDP metagenomic validation pipeline. 

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

#Check that input orthofile is present in the run directory
if [[ ! -f "$ORTHOFILE" ]]; then
    echo "The orthofinder file could not be opened:"
    echo $ORTHOFILE
    echo "Exiting script"
    exit 1
fi

#Check that mafft is installed and in path
if ! [ -x "$(command -v mafft)" ]; then
  echo 'Error: mafft is not installed/not in path.' >&2
  exit 1
fi

#Generate ortholog-file with single-copy ortholog families. Filter single-copy orthologs by length (gene-families containing genes shorter than 300bp are removed)
python3 bin/get_scp_orthologs.py $ORTHOFILE 
python3 bin/filter_orthologs.py single_ortho.txt

#Extract sequences of single-copy ortholog families into $outdir, move to $outdir
if [ -d $OUTDIR ]; then
    echo "Output directory already exists! Exiting script"
    exit 1
fi
echo "Extracting filtered single-copy ortholog gene-family sequences into directory: $OUTDIR"
python3 bin/extract_orthologs.py $OUTDIR
cd $OUTDIR

#For each gene-family: 1. align the amino-acid sequences, 2. translate to codon-aligned nucleotide alignement, 3. trim the alignment, 4. simplify seq-headers to genome-id
echo "Starting alignment of amino-acid sequences with mafft, back-translation and trimming"
for i in $( ls *.faa ); do
    OG=${i:0:9}
    echo "Processing: "$OG
    mafft --auto --quiet $OG".faa" > $OG"_aln.fasta"
    ALN_AA=$OG"_aln.fasta"
    FFN_FILE=$OG".ffn"
    python3 ../bin/aln_aa_to_dna.py $ALN_AA $FFN_FILE
    ALN_NUC=$OG"_aln_nuc.fasta"
    python3 ../bin/trim_aln.py $ALN_NUC
    ALN_TRIM=$OG"_aln_trim.fasta"
    sed 's/_.*$//g' $ALN_TRIM > temp_aln.txt
    mv temp_aln.txt $ALN_TRIM
done

#Generate a concatenate alignment in phylip-format, ready for phylogenetic inference with e.g. RAxML

python3 ../bin/concatenate_aln.py
echo "Generating concatenate alignment"
echo "All done"
