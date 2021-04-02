Bacterial species validation with metagenomic data
=======

This repository contains a pipeline for validating candidate bacterial species with metagenomic data.

If you are using the pipeline, please cite:

> Kirsten Maren Ellegaard & Philipp Engel. **Genomic diversity landscape of the honey bee gut microbiota**; _Nature Communications_ **10**, Article number: 446 (2019).
> PMID: 30683856;
> doi:[10.1038/s41467-019-08303-0](https://www.nature.com/articles/s41467-019-08303-0)

> Kirsten Maren Ellegaard, Shota Suenami, Ryo Miyasaki, Philipp Engel. **Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species**; _Current Biology_ **10**, Epub 2020 Jun 11.
> PMID: 32531278;
> doi: [10.1016/j.cub.2020.04.070](https://www.cell.com/current-biology/fulltext/S0960-9822(20)30586-8)
 
About metagenomic species validation: what and why
----------

The starting point of this pipeline is a collection of sequenced bacterial genomes belonging to the same 16S rRNA phylotype (i.e. > 97% identity in 16S rRNA gene). 16S rRNA phylotypes often contain multiple SDPs ("Sequence-discrete populations"), also  commonly referred to as "species", which by definition are discrete from each other at the sequence-level (see fx. [10.1038/s41467-018-07641-9](https://www.nature.com/articles/s41467-018-07641-9)). By analyzing bacterial genomes, one can get a good indication of whether a 16S rRNA phylotype contains multiple species. However, it is often unclear how well sequenced genomes represent natural bacterial populations, especially if only a small number of strains have been sequenced. If  genomes of bacterial isolates are to be used as a database for metagenomic analysis, it is therefore worthwhile to check whether:

* Putative related species in the database are discrete from each other in the metagenomes
* Genomes in the database are representative of the strains in the metagenomes 

The current pipeline will address both of these points.

The species validation is done based on core genes, i.e. genes known to be shared among all strains of the related species that are to be tested for discreteness. Core genes are widely used for inferring relatedness among bacteria, but also for estimating species abundances with metagenomic data. By using core genes for metagenomic species quantification, robust estimates may be obtained, even if the species engage in horizontal gene transfer for other genes, as is likely to be the case for closely related species that co-exist in the same environment. In the case of the honey bee gut microbiota, previous studies have used the largest possible set of core genes, i.e. the gene-set shared among all strains of a given 16S rRNA phylotype (as estimated from the available sequenced genomes). However, the pipeline can also be run with a smaller set of genes, i.e. universal core genes. This may be preferable if a large number of species are to be tested, or if there are very few available reference genomes. 

Pre-requisites
--------

This pipeline requires:

* Python 3 (or higher)
* numpy
* biopython
* Bash
* blast+ suite ([link]())
* muscle ([link]())


Installation
--------

```bash
git clone https://github.com/kirsten2/SDP_validation.git
```

Note: ```blast``` and ```muscle``` must be in the system path when executing the bash-scripts

Quick-start: Running the pipeline with example data
--------

The pipeline requires the following input-files:

1. Two directories (```faa_files```,```ffn_files```), containing gene-sequences for all reference genomes
2. A file with metagenomic ORFs predicted on metagenome assemblies
3. A tab-delimited file, specifying genome-ids for the candidate species
4. A file with orthologous gene-families, predicted on the reference genomes

A data-set derived from the honey bee gut microbiota can be downloaded from zenodo, and used to test the pipeline [link](https://sandbox.zenodo.org/record/710401#.X9h27i3Mx2c). Download the data-set from zenodo:

```bash
python3 bin/download_data.py --genome_db --metagenomic_orfs
```

**Expected result**: four directories with genomic data for all genomes in the honey bee gut microbiota database (```faa_files```,```ffn_files```,```bed_files```, ```gff_files```), the genome database (```genome_db_210402```), the genome database metafile (```genome_db_metafile_210402.txt```) and a fasta-file with metagenomic ORFs (```metagenomic_orfs.ffn```)

To run the pipeline on the 16S rRNA phylotype "Firm5", using 10 universal core gene families, start by generating alignment files using the reference genes of the sequenced genomes:

```bash
bash bin/prep_corefam_aln.sh -d firm5_uni -o OrthologousGroups_uni_example.txt
```

**Expected result**: a new directory named ```firm5_uni```, containing fasta-files for the sequences and alignments of the 10 universal core gene families. 

Next, run the species validation pipeline:

```bash
bash bin/species_validation.sh -c Candidate_species_uni_example.txt -i firm5_uni -d metagenomic_orfs.ffn
```

**Expected result**: Six new directories (```firm5_1```, ```firm5_4```, ```firm5_7```,), corresponding to each of the seven putative species affiliated with the 16S rRNA phylotype. Each directory contains fasta-files with sequences of the ORFs recruited to each core gene family, and a file named ```perc_id.txt``` with the alignment results. The file ```log.txt``` will be printed in the run-directory, containing some summary data on the results. 

The ```perc_id.txt```files contains the maximum alignment percentage identity for each recruited ORF to the reference core sequences of the recruiting species. It also details whether the first blast-hit for the ORF is to the recruiting species or to another species within the same 16S rRNA phylotype.

The R-script in the bin-directory can be used to plot the density distribution of the data, for example:

```bash
cd firm5_3
Rscript ../bin/plot_validation.R "firm5_3" perc_id.txt
```

**Expected result**: in this case, a file named ```firm5_3_recruitment_plot.pdf```within the ```firm5_3```species directory. 

To run the pipeline with the full set of core genes for all the putative species of this phylotype, re-run the bash-scripts substituting the candidate species file and the ortholog-file:


```bash
grep firm5 genome_db_metafile_210402.txt > Candidate_species_all_example.txt
bash bin/prep_corefam_aln.sh -d firm5_all -o Orthofinder/4_firm5_single_ortho_filt.txt
bash bin/species_validation.sh -c Candidate_species_all_example.txt -i firm5_all -d metagenomic_orfs.ffn
```

Check the wiki for further details on the pipeline and plot interpretation:

[wiki](https://github.com/kirsten2/SDP_validation/wiki) 

Running the pipeline with other data
---------

To run the pipeline with your own data, further instructions on file-preparation and formatting are provided in the wiki:

[wiki](https://github.com/kirsten2/SDP_validation/wiki) 