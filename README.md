This repository contains scripts for validating candidate SDPs with metagenomic data, using the approach described in the publications:

Ellegaard & Engel (2019): Genomic diversity landscape of the honey bee gut microbiota. Nature communications.
Ellegaard et al. (2020): Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species. Current Biology

The starting point of the pipeline is a collection of sequenced genomes (bacterial isolates), belonging to the same 16S rRNA phylotype (i.e. > 97% identity in 16S rRNA gene).

Briefly, the pipeline consists of the following steps:
1. Inference of single-copy core gene families
2. Construction of trimmed nucleotide alignments for each core gene family
3. Determination of candidate SDPs, based on core genome phylogenies and genomic average nucleotide identities
4. Validation of candidate SDPs, based on recruitment of core ortholog sequences from metagenomic ORFs
