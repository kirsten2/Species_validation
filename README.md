Bacterial species validation with metagenomic data
=======

This repository contains scripts for validating candidate bacterial species with metagenomic data.

If you are using this pipeline, please cite:

> Kirsten Maren Ellegaard & Philipp Engel. **Genomic diversity landscape of the honey bee > gut microbiota**; _Nature Communications_ **10**, Article number: 446 (2019).
> PMID: 30683856;
> doi:[10.1038/s41467-019-08303-0](https://www.nature.com/articles/s41467-019-08303-0)

> Kirsten Maren Ellegaard, Shota Suenami, Ryo Miyasaki, Philipp Engel. **Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species**; _Current Biology_ **10**, Epub 2020 Jun 11.
> PMID: 32531278;
> doi: [10.1016/j.cub.2020.04.070](https://www.cell.com/current-biology/fulltext/S0960-9822(20)30586-8)
 
About metagenomic species validation: what and why
----------

The starting point of this pipeline is a collection of sequenced bacterial genomes belonging to the same 16S rRNA phylotype (i.e. > 97% identity in 16S rRNA gene). 16S rRNA phylotypes often contain multiple SDPs ("Sequence-discrete populations"), also  commonly referred to as "species", which by definition are discrete from each other at the sequence-level (see fx. [10.1038/s41467-018-07641-9](https://www.nature.com/articles/s41467-018-07641-9)). 

By analyzing bacterial genomes, one can get a good indication of whether a 16S rRNA phylotype contains multiple species. However, it is often unclear how well sequenced genomes represent natural bacterial populations, especially if only a small number of strains have been sequenced. If  genomes of bacterial isolates are to be used as a database for metagenomic analysis, it is therefore worthwhile to check whether:

* Putative species in the database are discrete from each other in the metagenomes
* Genomes in the database are representative of the strains in the metagenomes 

The current pipeline will address both of these points.

The species validation is done based on the 16S rRNA phylotype core genes, i.e. the gene-set shared among all strains of the 16S rRNA phylotype, as inferred from the sequenced genomes. The rationale for using these genes for the validation is that they are also commonly used for estimating relative abundances of species with metagenomic data. If the species validation pipeline confirms that the species are discrete from each other based, on the core genes, we can also assume that they can be reliably quantified with these genes. This will be true, even if the species engage in horizontal gene transfer for other genes, as is likely to be the case for closely related species that co-exist in the same environment.

Pre-requisites
--------

This pipeline requires:

* Python 3 (or higher)
* numpy 1.15.0
* biopython..
* Bash (version 4 or higher)
* blast+ suite ([link]())
* muscle ([link]())


Installation
--------

```bash
git clone https://github.com/kirsten2/SDP_validation.git
```

Note: ```blast``` and ```muscle``` must be in the system path when executing the bash-scripts

Running pipeline with example data
--------

In this example, the pipeline is applied to three species belonging to the same 16S rRNA phylotype (> 97% 16S rRNA), using a genome database for the honey bee gut microbiota and two metagenomic samples ([zenodo_link](coming here)). 

Download example data from zenodo:

```bash
python bin/download_data.py --genes --metagenomes
```
**Expected result**: two directories with gene-sequences for all genomes in the honey bee gut microbiota database (```genes_faa```,```genes_ffn```), and a fasta-file with all metagenomic ORFs (```metagenomic_orfs.ffn```)

Subset gene-files for single-copy core gene family sequences (based on the example ortholog-prediction file), and generate alignment-files:

```bash
bash bin/prep_corefam_aln.sh -d firm5 -o OrthologousGroups_example.txt
```

**Expected result**: a new directory named ```firm5```, containing core gene family sequences and alignments. Estimated time: ~10min

Run the species validation pipeline:

```bash
bash bin/species_validation.sh -c Candidate_species_example.txt -i firm5 -d metagenomic_orfs.ffn
```

**Expected result**: Three new directories (```firm5_3```,```firm5_5```,```firm5_7```), containing fasta-files with sequences of ORFs recruited to each core gene family, and a file named ```perc_id.txt``` within each directory. The file ```log.txt``` will be printed in the run-directory, containing some summary data on the results. Estimated time: ~15min 

The ```perc_id.txt```files contains the maximum alignment percentage identity for each recruited ORF to the reference core sequences of the recruiting species. It also details whether the first blast-hit for the ORF is to the recruiting species or to another species within the same 16S rRNA phylotype.

The R-script in the bin-directory can be used to plot the density distribution of the data, for example:

```bash
cd firm5_3
Rscript ../bin/plot_validation.R "firm5_3" perc_id.txt
```

**Expected result**: in this case, a file named ```firm5_3_recruitment_plot.pdf```within the ```firm5_3```species directory. 

Check the wiki for further details on the pipeline and plot interpretation:

[wiki](https://github.com/kirsten2/SDP_validation/wiki) 

Running the pipeline with other data
---------

Further details on the pipeline are provided in the wiki:

[wiki](https://github.com/kirsten2/SDP_validation/wiki) 