NOTE: this repository is currently under development!

This repository contains scripts for validating candidate SDPs with metagenomic data, using the approach described in the publications:

Ellegaard & Engel (2019): Genomic diversity landscape of the honey bee gut microbiota. Nature communications.
Ellegaard et al. (2020): Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species. Current Biology

#About SDP validation: what and why

The starting point of the pipeline is a collection of sequenced genomes (bacterial isolates), belonging to the same 16S rRNA phylotype (i.e. > 97% identity in 16S rRNA gene). 16S rRNA phylotypes often contain multiple SDPs (aka "species"), which by definition are discrete from each other at the sequence-level (all strains of the same species are more closely related to each other than to strains of any other species). By analyzing genomes of bacterial isolates, one can get a good indication of whether a 16S rRNA phylotype contains multiple SDPs (see fx. Jain et al. 2019). However, working with bacterial isolates, it is usually unknown how well the sequenced genomes represent the natural bacterial population, especially if only a small number of strains have been sequenced. Therefore, I refer to such SDPs as "Candidate SDPs". 

If  genomes of bacterial isolates are to be used as a database for metagenomic analysis, it is therefore worthwhile to check whether the candidate SDPs are discrete from each other in the metagenomic samples. Moreover, if the metagenomic samples contain SDPs that lack representative sequenced bacterial isolates, then reads from such SDPs may interfere with the quantification of other SDPs. The current pipeline will address both of these points.

The SDP validation in this pipeline is done based on the 16S rRNA phylotype core genes (i.e. the gene-set shared among all strains of the 16S rRNA phylotype). A major motivation for using these genes is that they are also used for our metagenomic community profiling pipeline (repository in prep). Thus, if SDP validation pipeline confirms that the SDPs are discrete from each other based on the core genes, we can also be confident that we can quantify them reliably using these genes. This will be true, even if the SDPs engage in horizontal gene transfer in other genes, as is likely to be the case for closely related SDPs that co-exist in the same environment.

#Overview of pipeline steps

Briefly, the pipeline consists of the following steps:
1. Inference of single-copy core gene families (ortholog prediction), using isolate genomes
2. Construction of nucleotide alignments for each core gene family
3. Determination of candidate SDPs, based on core genome phylogenies and genomic average nucleotide identities
4. Validation of candidate SDPs, based on alignment of metagenomic ORFs to single-copy core gene families

The first two steps can be performed using the bash-script "corefam_aln.sh". It takes as input a file of orthologous gene-families, formatted as per Orthofinder. Furthermore, amino-acid and nucleotide sequences of the genes in the orthofinder file must be present in the directories "genome_faa/genome_ffn", which in turn must be present in the run directory.

The third step typically involves generation of core genome phylogenies and/or estimation of genomic ANI. For convenience. A concatenate alignment of core gene sequences is generated with the bash-script  "corefam_aln.sh", which can be used directly as input for many phylogenetic tools (fx. RAxML). However, I leave it up to the user to decide which approach to use here. At the end, the user will need to generate a plain text-file indicate the genome-ids of a set of candidate-SDPs.

The fourth step can be done with the bash-script "SDP_validation.sh". The overall idea is, as stated,  to align metagenomic ORFs to single-copy core gene families, and record the percentage identity for each ORF (last part still in prep.)

