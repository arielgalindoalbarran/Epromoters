# Epromoters


Here we present a collection of scripts used for the publication:

# Genome-wide characterization of mammalian promoters with distal enhancer functions.

Authors: Lan T.M. Dao, Ariel O Galindo-Albarrán, Jaime A Castro-Mondragon, Charlotte Andrieu-Soler, Alejandra Medina-Rivera, Charbel Souaid, Guillaume Charbonnier, Aurélien Griffon, Laurent Vanhille, Tharshana Stephen, Jaafar Alomairi, David Martin, Magali Torres, Nicolas Fernandez, Eric Soler, Jacques van Helden, Denis Puthier, Salvatore Spicuglia.

The scripts and the contribution to the results are listed:

1.	Human CapStarr-seq. 
	File: CapStarr-Seq.
	Include: Snakefile and code.

2.	RNA transcription and selection of control set.
	File: ChIP-Seq_and_DNAse-Seq and TF_RNA_Epigenomic.
	Include: R scripts, customized program coded in python: gtftoolkit.

3.	Epigenomic analysis, and Transcription factors enrichment and density.
	File: ChIP-Seq_and_DNAse-Seq and TF_RNA_Epigenomic.
	Include: R scripts, customized program coded in python: gtftoolkit.

4.	Motif analysis in Epromoters.
	File: Epromoters_motif_analysis.
	Include: R scripts and data. 
5.	Computations for ChIA-PET enrichment score of P-P interactions.
	File: ChIA-PET
	Include: R scripts, customized program coded in python: gtftoolkit 
6.	eQTL analysis.
	File: eQTL
	Include: R scripts, data and results.

Note: A [Conda](https://conda.io/docs/) installation with [Bioconda](https://bioconda.github.io/) channels is recommended to handle software dependencies.
The CapStarr-Seq, ChIA-PET, ChIP-Seq and DNAse-Seq analyses require to have [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed through Conda.


