#!/bin/bash

#Script purpose and preparation of file directory TO BE ADDED HERE

#Usage: bash bin/corefam_aln.sh -o OrthologousGroups.txt -d firm5

#Parse script arguments (name of ortholog-file and output directory).
while getopts o:d: flag
do
    case "${flag}" in
        o) orthofile=${OPTARG};;
        d) outdir=${OPTARG};;
    esac
done

#Check that input orthofile is present in the run directory
if [[ ! -f "$orthofile" ]]; then
    echo "$orthofile does not exist!"
    echo "Exiting script"
    exit
fi

#Check that mafft is installed and in path
if ! [ -x "$(command -v mafft)" ]; then
  echo 'Error: mafft is not installed/not in path.' >&2
  exit 1
fi

#Generate ortholog-file with single-copy ortholog families
python3 bin/get_scp_orthologs.py $orthofile 

#Extract sequences of single-copy ortholog families into $outdir, move to $outdir
echo "Extracting single-copy ortholog sequences into directory: $outdir"
python3 bin/extract_orthologs.py $outdir
cd $outdir

#For each gene-family: 1. align the amino-acid sequences, 2. translate to codon-aligned nucleotide alignement, 3. trim the alignment, 4. simplify seq-headers to genome-id
for i in $( ls *.faa ); do
    OG=${i:0:9}
    echo "Align amino-acid sequences with mafft, back-translate and trim: "$OG
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
echo "Script completed!"
